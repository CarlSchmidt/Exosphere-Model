;+
; NAME:
;	Irrandiance
;
; PURPOSE:
;	Provide the Solar spactral irradiance, time shifted and scaled 
;	to the solar distance of any major body from Earth-based measurements
;
; CALLING SEQUENCE:  
;	Irrandiance, body, UTC, Flux
;
; INPUTS:
;	body			  Body Name String 
;	UTC			    Universal Coordiate Time String
;
; OUTPUTS:  
; Flux        2xN array for N wavelength bins, wavelength in nm, flux in photons cm^-2 a^-1 
;
; COMMON BLOCKS:
;	
; PROCEDURE:
;	1.  Check input parameters, Set up error handing
;	2.  Calculate position of Earth and body in heliocentric inertial frame
;	3.  Calculate solar rotation angle between Earth and given body
;	    and determine optimal number of days for time shifting 
;	    Earth-based irradiance to the body
;	4.  Pull the TIMED / SEE Flux and the SORCE flux from the LISIRD API for the shifted times
;	    Apply weighting of 2 Carrington rotations if more than 1/4 of a solar rotation seperates the measurement times.
;	    

; MODIFICATION HISTORY:
;	Carl Schmidt, Univ. of Virginia, 2013

;+

pro Irrandiance, body, UTC, Time_range, directory, Flux, debug=debug

;	1.  Check input parameters 

  if (n_params() lt 3) then begin
    print, 'USAGE:  Irrandiance, body, UTC, Flux
    print, ' '
    print, '  body        Input Body Name String 
    print, '  UTC         Input Universal Coordinate Time String
    print, ' '
    print, '  Flux        Output 2xN array for N wavelength bins, wavelength in nm, flux in photons cm^-2 a^-1 
    print, ' '
  endif

  CATCH, Error_status
  IF Error_status NE 0 THEN BEGIN
      PRINT, 'Irradiance.pro Error index: ', Error_status
      PRINT, 'Irradiance.pro Error message: ', !ERROR_STATE.MSG
      print, '--> Using the default date of Jan 1, 2010 for solar drivers' 
      restore, strcompress(directory + '\Solar_Spectral_Irradiance\Default_Solar_Flux.sav')
      goto, use_default_flux
      CATCH, /CANCEL
  ENDIF
  
  ;If the irradiance was calculated for this body, at this time, in the last run then skip to the end to save time
  restore, strcompress(directory + 'Solar_Irradiance.sav')
  if ((Solar_Irradiance.UTC eq UTC) and (Solar_Irradiance.body eq body)) then begin
    print, strcompress('Loading saved flux at ' + Solar_Irradiance.body + ' on ' + Solar_Irradiance.UTC +' . . .')  
    flux = Solar_Irradiance.flux
    goto, load_flux   
  endif
  
;	2.  Calculate heliocentric position of Earth and planet:
  
  cspice_str2et, UTC, ephemeris_time; convert UTC to ephemeris time (expressed as the number of ephemeris seconds past J2000)  
  
  ;Get the coordinates in the Helicentric Inertial frame, that is:
    ;
    ;    -  X-Y plane is defined by the Sun's equator of epoch J2000: the +Z
    ;       axis, primary vector, is parallel to the Sun's rotation axis of
    ;       epoch J2000, pointing toward the Sun's north pole;
    ;
    ;    -  +X axis is defined by the ascending node of the Sun's equatorial
    ;       plane on the ecliptic plane of J2000;
    ;
    ;    -  +Y completes the right-handed frame;
    ;
    ;    -  the origin of this frame is the Sun's center of mass.
    
    ;Since the transformation between the HCI frame and J2000 frame is fixed and time independent, converting between them
    ;doesn't need a specialized SPICE kernel. See NAIF .TF files in the STEREO frames kernel for more info:
   
  J2000_to_HCI_MATRIX        = [$
            [ 0.2458856764679510,       0.8893142951159845,       0.3855649343628876],$
            [-0.9615455562494245,       0.1735802308455697,       0.2128380762847277],$
            [ 0.1223534934723278,      -0.4230720836476433,       0.8977971010607901]]
  
  ;EARTH
    cspice_spkpos, 'Earth' , ephemeris_time, 'J2000', 'LT', 'Sun', state, light_time ;state is [x,y,z] with respect to the Sun
    state = transpose(J2000_to_HCI_MATRIX) # state  
    cspice_reclat, state, Earth_distance, Earth_longitude, Earth_latitude
    
  ;BODY 
    cspice_spkpos, Body, ephemeris_time, 'J2000', 'LT', 'Sun', state, light_time ;state is [x,y,z] with respect to the Sun
    state = transpose(J2000_to_HCI_MATRIX) # state  
    cspice_reclat, state, Body_distance, Body_longitude, Body_latitude 

;	3.  Calculate solar rotation angle between Earth and given planet and determine optimal number of days for time shifting 
;	    Earth-based irradiance to the planet

  angle = (Earth_longitude - Body_longitude)*!RADEG   ; in degrees
  if (angle gt 180.) then angle = angle - 360.
  if (angle lt -180.) then angle = angle + 360.
  days_rotate = 27.0                                  ; Assumes 27-days for solar rotation (360 degrees)
  flt_shift =  days_rotate * angle / 360.
  
  ;	Make time_shift_days (tsd) a 2 x 2 array
  ;	     time_shift_days[*,0] = + / - shift days [ -(days ago) , (days ahead)] 
  ;      time_shift_days[*,1] = weights in making average
  ;
  ;	Weight 7 - 13.5 day sfits, but force a weighting of 1.0 for days less than 7 days (i.e., trust values within 1/4 solar rotation)

  shift1 = flt_shift
  if (shift1 lt 0) then shift2 = days_rotate + shift1 else shift2 = shift1 - days_rotate
  wt1 = (0.5 + (days_rotate/2. - abs(shift1)) * 0.077) < 1.0
  wt2 = 1. - wt1
  time_shift_days = [ [shift1, shift2], [wt1, wt2] ]

;4: Pull the TIMED / SEE Flux and the SORCE flux from the LISIRD API for the shifted times, use the time at the mid-point of the integration 
    
    ;Dates prior to 2002 or after 2013 don't have solar flux measurements; use a default value in this case
    cspice_str2et, SYSTIME(/UTC), current_time
    if (ephemeris_time lt 1.0242726e+008) or (ephemeris_time gt current_time) then begin 
      print, 'No solar flux measurements exist for the date ', UTC
      print, '--> Using the default date of Jan 1, 2010 for solar drivers' 
      restore, strcompress(directory + '\Solar_Spectral_Irradiance\Default_Solar_Flux.sav')
      goto, use_default_flux
    endif
    
    maxtime = max(Time_range)*24.*3600. ;change model's stop time from days to seconds
    mintime = min(Time_range)*24.*3600. ;change model's start time from days to seconds
    
    no_weighting = where(time_shift_days[*,1] eq 1)
    if (no_weighting ne -1) then begin
        shift = time_shift_days[no_weighting,0] * 24.*3600.

        ;get the day of year format for this time.
        cspice_et2utc, ephemeris_time + shift - ((maxtime-mintime)/2. + mintime), 'ISOC', 0, Irradiance_date  
        
        ;get the irrandiance from SEE using the LASP LISIRD API
        SEE = webget(strcompress('http://lasp.colorado.edu/lisird/tss/timed_see_ssi.csv?wavelength,irradiance&time~' + strmid(Irradiance_date, 0, 10)))
          result = STRSPLIT(STRJOIN(SEE.text, ' ', /single), ' ', /extract)
          openw, lun, strcompress(directory + 'Solar_Spectral_Irradiance\SEE.txt'), /get_lun
          for i = 4, n_elements(result)-2 do begin
            printf, lun, result[i]
          endfor
          close, lun
          free_lun, lun
        READCOL,strcompress(directory + 'Solar_Spectral_Irradiance\SEE.txt'), Format='A,A', wavelength, intensity, /SILENT  
        start = where(wavelength eq '0.5')
        SEE_wavelength = double(wavelength[start:*])                    ;in nm 
        SEE_intensity = double(intensity[start:*]) * 1.e3               ;flux in mW/m^2/nm (using a W to mW conversion)
  
        ;get the irrandiance from SORCE using the LASP LISIRD API
          SORCE = webget(strcompress('http://lasp.colorado.edu/lisird/tss/sorce_ssi.csv?wavelength,irradiance&time~' + strmid(Irradiance_date, 0, 10)))
          ;Write the SORCE data to an ASCII text file
            result = STRSPLIT(STRJOIN(SORCE.text, ' ', /single), ' ', /extract)
            openw, lun, strcompress(directory + 'Solar_Spectral_Irradiance\SORCE.txt'), /get_lun
            for i = 4, n_elements(result)-2 do begin
              printf, lun, result[i]
            endfor
            close, lun
            free_lun, lun
          READCOL,strcompress(directory + 'Solar_Spectral_Irradiance\SORCE.txt'), Format='A,A', wavelength, intensity, /SILENT  
          start = where(wavelength eq '0.50')
          SORCE_wavelength = double(wavelength[start:*])                    ;in nm 
          SORCE_intensity = double(intensity[start:*]) * 1.e3               ;flux in mW/m^2/nm (using a W to mW conversion)
          If keyword_set(debug) then begin
            plot, SORCE_wavelength, SORCE_intensity, xrange = [0, max(SORCE_wavelength)], /ylog, yrange = [1.e-5, 1.e4], ystyle = 1., Xtitle = 'Wavelength (nm)', $
                ytitle = 'Solar Irradiance (W/m^2/nm)', Title = 'Incident Solar Flux'
            oplot, SEE_wavelength, SEE_intensity, linestyle = 2. ;overlap looks okay? 
          endif

        sorce = [transpose(SORCE_wavelength[where(SORCE_wavelength eq 190.5):*]), transpose(SORCE_intensity[where(SORCE_wavelength eq 190.5):*])]
        see = [transpose(see_wavelength[0:where(see_wavelength eq '189.50')]),transpose(SEE_intensity[0:where(see_wavelength eq '189.50')])]
        flux = [[see], [sorce]]
    endif else begin ;weighted average of two solar rotations
      shifts = time_shift_days[*,0] * 24.*3600.
      
      ;What was the solar flux was at the previous pass?
        ;get the day of year format for this time.
        cspice_et2utc, ephemeris_time + shifts[0] - ((maxtime-mintime)/2. + mintime), 'ISOC', 0, Irradiance_date  
        
        ;get the irrandiance from SEE using the LASP LISIRD API
        SEE = webget(strcompress('http://lasp.colorado.edu/lisird/tss/timed_see_ssi.csv?wavelength,irradiance&time~' + strmid(Irradiance_date, 0, 10)))
          result = STRSPLIT(STRJOIN(SEE.text, ' ', /single), ' ', /extract)
          openw, lun, strcompress(directory + 'Solar_Spectral_Irradiance\SEE_before.txt'), /get_lun
          for i = 4, n_elements(result)-2 do begin
            printf, lun, result[i]
          endfor
          close, lun
          free_lun, lun
        READCOL, strcompress(directory + 'Solar_Spectral_Irradiance\SEE_before.txt'), Format='A,A', wavelength, intensity, /SILENT  
        start = where(wavelength eq '0.5')
        SEE_before_wavelength = double(wavelength[start:*])                    ;in nm 
        SEE_before_intensity = double(intensity[start:*]) * 1.e3               ;flux in mW/m^2/nm (using a W to mW conversion)
  
        ;get the irrandiance from SORCE using the LASP LISIRD API
          SORCE = webget(strcompress('http://lasp.colorado.edu/lisird/tss/sorce_ssi.csv?wavelength,irradiance&time~' + strmid(Irradiance_date, 0, 10)))
          ;Write the SORCE data to an ASCII text file
            result = STRSPLIT(STRJOIN(SORCE.text, ' ', /single), ' ', /extract)
            openw, lun, strcompress(directory + 'Solar_Spectral_Irradiance\SORCE_before.txt'), /get_lun
            for i = 4, n_elements(result)-2 do begin
              printf, lun, result[i]
            endfor
            close, lun
            free_lun, lun
          READCOL,strcompress(directory + 'Solar_Spectral_Irradiance\SORCE_before.txt'), Format='A,A', wavelength, intensity, /SILENT  
          start = where(wavelength eq '0.50')
          SORCE_before_wavelength = double(wavelength[start:*])                    ;in nm 
          SORCE_before_intensity = double(intensity[start:*]) * 1.e3               ;flux in mW/m^2/nm (using a W to mW conversion)
      
      ;What will the flux be on the next pass? 
        ;get the day of year format for this time.
        cspice_et2utc, ephemeris_time + shifts[1] - ((maxtime-mintime)/2. + mintime), 'ISOC', 0, Irradiance_date  
        
        ;get the irrandiance from SEE using the LASP LISIRD API
        SEE = webget(strcompress('http://lasp.colorado.edu/lisird/tss/timed_see_ssi.csv?wavelength,irradiance&time~' + strmid(Irradiance_date, 0, 10)))
          result = STRSPLIT(STRJOIN(SEE.text, ' ', /single), ' ', /extract)
          openw, lun, strcompress(directory + 'Solar_Spectral_Irradiance\SEE_after.txt'), /get_lun
          for i = 4, n_elements(result)-2 do begin
            printf, lun, result[i]
          endfor
          close, lun
          free_lun, lun
        READCOL, strcompress(directory + 'Solar_Spectral_Irradiance\SEE_after.txt'), Format='A,A', wavelength, intensity, /SILENT  
        start = where(wavelength eq '0.5')
        SEE_after_wavelength = double(wavelength[start:*])                    ;in nm 
        SEE_after_intensity = double(intensity[start:*]) * 1.e3               ;flux in mW/m^2/nm (using a W to mW conversion)
  
        ;get the irrandiance from SORCE using the LASP LISIRD API
          SORCE = webget(strcompress('http://lasp.colorado.edu/lisird/tss/sorce_ssi.csv?wavelength,irradiance&time~' + strmid(Irradiance_date, 0, 10)))
          ;Write the SORCE data to an ASCII text file
            result = STRSPLIT(STRJOIN(SORCE.text, ' ', /single), ' ', /extract)
            openw, lun, strcompress(directory + 'Solar_Spectral_Irradiance\SORCE_after.txt'), /get_lun
            for i = 4, n_elements(result)-2 do begin
              printf, lun, result[i]
            endfor
            close, lun
            free_lun, lun
          READCOL,strcompress(directory + 'Solar_Spectral_Irradiance\SORCE_after.txt'), Format='A,A', wavelength, intensity, /SILENT  
          start = where(wavelength eq '0.50')
          SORCE_after_wavelength = double(wavelength[start:*])                    ;in nm 
          SORCE_after_intensity = double(intensity[start:*]) * 1.e3         ;flux in mW/m^2/nm (using a W to mW conversion)

      ;Before taking the weighted average of the two solar rotations, check to be sure that no bogus values exist in either instance.
      for i = 0, 3 do begin
         case i of
            0: flux = [transpose(SEE_before_wavelength), transpose(SEE_before_intensity)]
            1: flux = [transpose(SEE_after_wavelength), transpose(SEE_after_intensity)]
            2: flux = [transpose(SORCE_before_wavelength), transpose(SORCE_before_intensity)]
            3: flux = [transpose(SORCE_after_wavelength), transpose(SORCE_after_intensity)]
         endcase
         ;spikes in wavelength indicate something bad, remove if spike are present 
         if not array_equal(flux[0,*], flux(0,sort(flux[0,*]))) then begin
            result = max(abs(flux[0,*] - smooth(flux[0,*],3)), badwavelength)        ;bad wavelength = wavelngth array outlier 
            keep = cgSetDifference(indgen(n_elements(flux[0,*])), badwavelength)     ;indicies of good wavelengths
            flux = transpose([[transpose(flux[0,keep])], [transpose(flux[1,keep])]]) ;new shortened array with bad wavelengths tossed out
         endif
         bogus_wavelength = where((flux[0,*] eq 1) or (flux[0,*] eq 0) or (finite(flux[0,*]) eq 0), COMPLEMENT=goodwavelength, NCOMPLEMENT=goodcount) 
         bogus_intensity = where((flux[1,*] eq 1) or (flux[1,*] le 0.) or (finite(flux[1,*]) eq 0), COMPLEMENT=goodintensity, NCOMPLEMENT=goodcount) 
         keep = cgSetIntersection(goodwavelength, goodintensity)
         flux = transpose([[transpose(flux[0,keep])], [transpose(flux[1,keep])]])
         case i of
            0: SEE_before_wavelength =   flux[0,*] 
            1: SEE_after_wavelength =    flux[0,*] 
            2: SORCE_before_wavelength = flux[0,*] 
            3: SORCE_after_wavelength =  flux[0,*] 
         endcase
         case i of
            0: SEE_before_intensity =   flux[1,*]
            1: SEE_after_intensity =    flux[1,*]
            2: SORCE_before_intensity = flux[1,*]
            3: SORCE_after_intensity =  flux[1,*]
         endcase
      endfor

    ;Combine a weighted average of the two Carrington rotations
    
      if array_equal(SEE_before_wavelength, SEE_after_wavelength) then begin
         ;Take the weighted average over all wavelengths
         SEE_wavelength = SEE_before_wavelength 
         SEE_intensity = (SEE_before_intensity*time_shift_days[0,1] + SEE_after_intensity*time_shift_days[1,1]) 
      endif else begin
         ;Take the weighted average, but only if a wavelength bin wasn't rejected (still exists) in both arrays 
         ;Find the wavelength bins that overlap in each (after some were rejected as bogus)  
         match, REFORM(SEE_before_wavelength), REFORM(SEE_after_wavelength), sub_SEE_before_wavelength, sub_SEE_after_wavelength, Count = Count_SEE_overlap
         if n_elements(SEE_before_wavelength) gt n_elements(SEE_after_wavelength) then begin
             SEE_wavelength = SEE_before_wavelength
             SEE_intensity = SEE_before_intensity
             SEE_intensity[sub_SEE_before_wavelength] = SEE_before_intensity[sub_SEE_before_wavelength]*time_shift_days[0,1] + $
                                                        SEE_after_intensity[sub_SEE_after_wavelength]*time_shift_days[1,1]
         endif else begin
             SEE_wavelength = SEE_after_wavelength
             SEE_intensity = SEE_after_intensity
             SEE_intensity[sub_SEE_after_wavelength] = SEE_before_intensity[sub_SEE_before_wavelength]*time_shift_days[0,1] + $
                                                       SEE_after_intensity[sub_SEE_after_wavelength]*time_shift_days[1,1]
         endelse
      endelse  
      if array_equal(SORCE_before_wavelength, SORCE_after_wavelength) then begin
         ;Take the weighted average over all wavelengths
         SORCE_wavelength = SORCE_before_wavelength
         SORCE_intensity = (SORCE_before_intensity*time_shift_days[0,1] + SORCE_after_intensity*time_shift_days[1,1]) 
      endif else begin
         ;Take the weighted average, but only if a wavelength bin wasn't rejected (still exists) in both arrays 
         ;Find the wavelength bins that overlap in each (after some were rejected as bogus) 
         match, REFORM(SORCE_before_wavelength), REFORM(SORCE_after_wavelength), sub_SORCE_before_wavelength, sub_SORCE_after_wavelength, Count = Count_SORCE_overlap        
         if n_elements(SORCE_before_wavelength) gt n_elements(SORCE_after_wavelength) then begin
             SORCE_wavelength = SORCE_before_wavelength
             SORCE_intensity = SORCE_before_intensity
             SORCE_intensity[sub_SORCE_before_wavelength] = SORCE_before_intensity[sub_SORCE_before_wavelength]*time_shift_days[0,1] + $
                                                            SORCE_after_intensity[sub_SORCE_after_wavelength]*time_shift_days[1,1]
         endif else begin
             SORCE_wavelength = SORCE_after_wavelength
             SORCE_intensity = SORCE_after_intensity
             SORCE_intensity[sub_SORCE_after_wavelength] = SORCE_before_intensity[sub_SORCE_before_wavelength]*time_shift_days[0,1] + $
                                                           SORCE_after_intensity[sub_SORCE_after_wavelength]*time_shift_days[1,1]
         endelse 
      endelse  

      ;take the weighted average of only the kept values
      see =   [transpose(SEE_wavelength[0:where(SEE_wavelength eq '189.50')]), $
               transpose(SEE_intensity[0:where(SEE_wavelength eq '189.50')])]
      sorce = [transpose(SORCE_wavelength[where(SORCE_wavelength eq 190.5):*]), $
               transpose(SORCE_intensity[where(SORCE_wavelength eq 190.5):*])]
      flux = [[see], [sorce]]
   endelse   
   use_default_flux: ;if no flux measurements exist, the program just skipped to this line and loaded a default solar flux on Jan 1, 2010 
      
      ;remove bogus values
         ;spikes in wavelength indicate something bad, remove if spike are present 
         if not array_equal(flux[0,*], flux(0,sort(flux[0,*]))) then begin
            result = max(abs(flux[0,*] - smooth(flux[0,*],3)), badwavelength)        ;bad wavelength = wavelngth array outlier 
            keep = cgSetDifference(indgen(n_elements(flux[0,*])), badwavelength)     ;indicies of good wavelengths
            flux = transpose([[transpose(flux[0,keep])], [transpose(flux[1,keep])]]) ;new shortened array with bad wavelengths tossed out
         endif
         bogus_wavelength = where((flux[0,*] eq 1) or (flux[0,*] eq 0) or (finite(flux[0,*]) eq 0), COMPLEMENT=goodwavelength, NCOMPLEMENT=goodcount) 
         bogus_intensity = where((flux[1,*] eq 1) or (flux[1,*] le 0.) or (finite(flux[1,*]) eq 0), COMPLEMENT=goodintensity, NCOMPLEMENT=goodcount) 
         keep = cgSetIntersection(goodwavelength, goodintensity)
         flux = transpose([[transpose(flux[0,keep])], [transpose(flux[1,keep])]])
      
      ;if the coverage is incomplete, fill in any gaps with the default solar flux on Jan 1, 2010.
      flux_on_UTC = flux
      restore, strcompress(directory + '\Solar_Spectral_Irradiance\Default_Solar_Flux.sav')
      default_flux = flux 
      match, reform(default_flux[0,*]), reform(flux_on_UTC[0,*]), sub_default_flux, sub_flux_on_UTC, Count = Count_overlap
      if count_overlap lt n_elements(default_flux[0,*]) then begin       
        default_flux[1,sub_default_flux] = flux_on_UTC[1,sub_flux_on_UTC]
        flux = default_flux
        print, strcompress('Solar flux on ' + UTC + ' only available for' + string(min(flux_on_UTC[0,sub_flux_on_UTC]), format = '(F6.1)') + ' to ' + $
          string(max(flux_on_UTC[0,sub_flux_on_UTC]), format = '(F6.1)') + ' nm. Using Jan 1, 2010 default flux for ' + $
          string(n_elements(default_flux[0,*]) - Count_overlap, format = '(I)') + ' other wavelength bins.') 
      endif else begin
        flux = flux_on_UTC 
      endelse  

      ;change flux units from mW/m^2/nm to photons/(cm^2 s nm)
      conversion = (1./1000.)*(1./100.^2.)/((6.62606957e-34*299792458.)/(flux[0,*]*1.e-9))
      flux[1,*] = flux[1,*] * conversion
      
      ;change flux units from mW/m^2/nm to ergs/(cm^2 s nm), They are equal!
      ;      conversion = (1./1000.) * (1.e7) * (1./100.^2.)
      ;      flux[1,*] = flux[1,*] * conversion
      
      ;scale flux with heliocentric distance
      flux[1,*] = flux[1,*] * (Earth_distance/Body_distance)^2.
      debug = 1
      If keyword_set(debug) then begin 
        print, ' '
        print, 'INPUT:  ', body, '  Time = ', UTC
        print, ' '
        print, 'Earth at (AU, Lon (deg), Lat (deg)', Earth_distance/149597871., Earth_longitude*!RADEG, Earth_latitude*!RADEG     
        print, 'Body at (AU, Lon (deg), Lat (deg) ', Body_distance/149597871., Body_longitude*!RADEG, Body_latitude*!RADEG
        print, 'Time Shift Used for Solar Flux (Days) & (Weights) = ', time_shift_days, format='(A,F8.3,"  [",F7.3,"]  &",F8.3,"  [",F6.3,"]")'
        Window, 0, Title = 'Incident Solar Flux',xs=600, ys=400 
        plot, flux[0,*]*10., flux[1,*]/10., xrange = [0, 24000.], /xstyle, Xtitle = strcompress('Wavelength (' + cgSymbol("angstrom")+')'), $
          ytitle = strcompress('Erg s!U-1!N cm!U-2!N ' + cgSymbol("angstrom") + '!U-1!N'), psym=10, $
          Title = strcompress('Solar Flux incident at ' + Body + ' at ' + UTC), charsize=1.6  
      endif 
      
      ;Save the output for a fast restore if this time and body are requested again.
      OpenW, lun, strcompress(directory + 'Solar_Irradiance_Photons.txt'), /get_lun
      PrintF, lun, flux, FORMAT='(F10.2,1X,E20.10)'
      Free_LUN, lun
      
      load_flux:
return
end

 
;====================================== NetCDF Alternative to TIMED / SEE API=======================================================

;;get the flux from SEE, use the time at the mid-piont of the integration  
;  cspice_et2utc, ephemeris_time - ((maxtime-mintime)/2. + mintime), 'ISOC', 0, Irradiance_date ;get the day of year format for this time. 
;  plot_see, round(date_conv(Irradiance_date)), planet = 'Earth', file = strcompress(directory + '\Solar_Spectral_Irradiance\' + Irradiance_data), $
;    type = 'SP', data = SEE ;flux is a daily average, rounded to the nearest day 
;  !psym = 0 ;LASP plot_see code sets the plot symbol variable, reset to default  
