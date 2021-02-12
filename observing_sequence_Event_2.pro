Pro Observing_Sequence_Event_2, Meteor_impact_UTC = Meteor_impact_UTC, Plume_Temperature = Plume_Temperature, $
                          Surface_distribution = Surface_distribution, loop_times = loop_times

; Runtimes:

; Background: 1.e26 1200K MBF, 100 loops, 0.166 Days duration = 2 days runtime (isotropic dayside?)

; Current Run
  Meteor_impact_UTC             = '2013-04-13 13:33:00'       ; time of the impact
  Plume_Temperature             = '10000K'                    ; temperature of the impact vapour
  Surface_distribution = 'Point_[235, 40]'                    ; Location of the impactor
  loop_times                    = 4.                          ; Bear minimum for any reasonable S/N (was 90 in Event 1)
  Na_Lofted                     = 1.e25                       ; seems like a lot
  Brightness_multiplier_Na      = 12.5   

COMMON Output_shared, Plot_range, Output_Size_In_Pixels, Output_Title, Center_in_frame, viewpoint, FOV, N_ticks, Tickstep, Observatory, Above_Ecliptic, Boresight_Pixel, Aperture_Corners
COMMON Model_shared, Body, Ephemeris_time, Seed, Directory, Particle_data, Line_data, Debug

Na_Data_Color  = 'Orange'
Na_Model_Color = 'Orange'
key_frames_Na  = [indgen(110)] * 8
;key_frames_Na  = [ [indgen(32)] * 8, [indgen(12)] * 40 + 256]
;key_frames_Na  = [indgen(51)] * 16
;key_frames_Na  = [indgen(27)] * 8
key_frames_Na  = key_frames_Na[UNIQ(key_frames_Na, SORT(key_frames_Na))]

  Na_scale = alog10([10^(-0.25), 10^1.75]) ; range to logrithmically scale the Na images in KiloRayleighs
  Mg_scale = alog10([10^(-1.5), 10^0.5]) ; range to logrithmically scale the Mg images in KiloRayleighs

; Calculate the flatness coefficient and planetary radius
  cspice_bodvrd, body, 'RADII', 3, radii
  re   = radii[0]
  rp   = radii[2]
  flat = ( re - rp ) / re

; Print the longitude and latitude of equatorial dawn at the impact time
  cspice_utc2et, Meteor_impact_UTC, meteor_ET
  cspice_subslr, 'Near point: ellipsoid', body, meteor_ET, 'IAU_'+body, 'none', viewpoint, sub_solar_point_planet_frame, trgepc, srfvec ; not sure why we'd need light time
  cspice_recpgr, body, sub_solar_point_planet_frame, re, flat, sub_solar_lon, sub_solar_lat, sub_solar_radius ; Convert to planetographic latitude and longitude, in the 'IAU_' body fixed frame  
  print, 'planetographic coordinates of equatorial dawn [Lon, Lat]', (sub_solar_lon + !pi/2. + 2.*!pi) mod (2.*!pi)*!radeg, sub_solar_lat*!radeg

; Load UVVS DDR from https://pds-geosciences.wustl.edu/messenger/mess-e_v_h-mascs-3-virs-cdr-caldata-v1/messmas_2101/data/ddr/ob2/uvvs_atmosphere/07/
  UVVS_DDR_Na      = read_mascs_ddr(Directory+'MESSENGER_UVVS\ud_09_ns_na.dat') 
  UVVS_UTC_TIME_Na = string(UVVS_DDR_Na.UTC_TIME)  ; convert time from byte to string YYDOYTHH:MM:SS.00

; define the time window to be plotted  
  ;plot_times       = ['2013-04-13 13:40', '2013-04-13 15:10'] ; times to plot (UTC)
  plot_times       = ['2013-04-13 13:33', '2013-04-13 15:10'] ; times to plot (UTC)
  cspice_utc2et, plot_times[0], start_ET
  cspice_utc2et, plot_times[1], stop_ET

; For Na: restrict to these times, format the dates in ephemeris time and Julian date
  SNR_threshold = 0.
  date      = UVVS_UTC_TIME_Na
  PDS_ET_Na = dblarr(N_elements(date))
  PDS_JD_Na = dblarr(N_elements(date))
  year      = '20' + StrMid(StrTrim(date,2), 0, 2)
  dayofyear = StrMid(StrTrim(date,2), 2, 3)
  time      = StrMid(StrTrim(date,2), 6, 11)
  CALDAT, JULDAY(1, dayofyear, year), month, day
  Na_UTC_string = year + '-' + string(month ,format = "(I2.2)") + '-' + string(day ,format = "(I2.2)") + ' ' + time
  for i = 0, N_elements(date)-1 do begin
    cspice_utc2et, Na_UTC_string[i], ET
    PDS_ET_Na[i] = ET
    cspice_et2utc, et, 'J', 6, JD_string
    PDS_JD_Na[i] = strmid(JD_string, 3)
  endfor
  duration      = where((PDS_ET_Na gt start_ET) and (PDS_ET_Na lt stop_ET), /NULL)
  keep_Na       = where(abs(UVVS_DDR_Na[duration].TOTAL_RADIANCE_SNR) gt SNR_threshold, /Null)
  DDR_Na_Subset = UVVS_DDR_Na[duration[Keep_Na]]
  DDR_Na        = UVVS_DDR_Na[duration[Keep_Na]].TOTAL_RADIANCE_KR
  DDR_Na_err    = UVVS_DDR_Na[duration[Keep_Na]].TOTAL_RADIANCE_KR / UVVS_DDR_Na[duration[Keep_Na]].TOTAL_RADIANCE_SNR
  PDS_ET_Na     = PDS_ET_Na[duration[Keep_Na]]
  PDS_JD_Na     = PDS_JD_Na[duration[Keep_Na]]
  Na_UTC_string = Na_UTC_string[duration[Keep_Na]]

; Load the background model's from Tim Cassidy and Matt Burger
  readcol, Directory+'MESSENGER_UVVS\Background_Model\Na_time_series_1884.txt', Na_BG_UTC_string_1884, Na_Obs_BG_1884, Na_BG_1884, format = 'A,F,F'
  Na_BG_UTC_string = Na_BG_UTC_string_1884
  Na_BG_JD = dblarr(N_elements(Na_BG_UTC_string))
  for i = 0, N_elements(Na_BG_UTC_string)-1 do begin
    cspice_utc2et, Na_BG_UTC_string[i], ET
    cspice_et2utc, et, 'J', 6, JD_string
    Na_BG_JD[i] = strmid(JD_string, 3)
  endfor  
  ; The times from these save files might be exactly the same as the PDS, so interpolate things to the PDS times
    Na_BG = interpol(Na_BG_1884, Na_BG_JD, PDS_JD_Na) 
    Na_BG[where(Na_BG le 0., /null)] = !Values.F_Nan

; Sometimes we want to simulate a metoer impact that occcured after the data sequence begins  
  key_frames_Na = key_frames_Na[where(PDS_ET_Na[key_frames_Na] gt meteor_ET, /NULL)]

; Format these inputs into a string and make a new write directory for this observing sequence  
  cspice_et2utc, meteor_ET, 'ISOC', 0, impact_UTC
  Surface_distribution_string = strmid(Surface_distribution, 7)
  Surface_distribution_string = strcompress(repstr(Surface_distribution_string,']','lat'), /remove_all)
  impact_UTC                  = repstr(impact_UTC,':','-') ; need to replace the colons in the directory string, use '-'  
  cd, directory
  write_directory             = impact_UTC + '---' + Surface_distribution_string + '---' + Plume_Temperature
  FILE_MKDIR, write_directory, /NOEXPAND_PATH

; Run the Na simulations at each observation time
  days_since_meteor_impact    = (PDS_ET_Na - meteor_ET) / (24.*3600.)

  for i = 0, N_elements(DDR_Na)-1 do begin
    junk = where(key_frames_Na eq i, count, /NULL)
    if count eq 0 then continue
    if days_since_meteor_impact[i] lt 0. then continue ; Can't model an event that hasn't happened yet    
    if FILE_TEST('C:\IDL\Generic Model V2\read_write\'+write_directory+'\Na_Meteor_frame_'+strcompress(string(i), /remove_all)+'.fit') then continue ; skip if we've done this before
    
    if i gt key_frames_Na[0] then restore_aloft_filename = 'C:\IDL\Generic Model V2\read_write\' + write_directory + '\Na_Meteor_frame_'+strcompress(key_frames_Na[junk-1], /remove)+'_loc_aloft.sav' else restore_aloft_filename = !Null
    
    Generic_Model, Time_range_this_run = days_since_meteor_impact[i], Output_title_this_run = write_directory + '\Na_Meteor_frame_'+strcompress(string(i), /remove_all), $
                   test_particle_this_run = 'Na', Line_this_run = 'Na-D', UTC_this_run = Na_UTC_string[i], $
                   Surface_distribution_this_run = Surface_distribution, Speed_distribution_this_run = 'MBF_'+plume_temperature, $
                   loop_times_this_run = loop_times, Upward_flux_at_exobase_this_run = Na_Lofted, restore_aloft_filename = restore_aloft_filename                     
  endfor

;-----MOVIE-----------------------------------------
  mpgFilename = Directory+Write_directory+'\Correlation_Sequence.mp4'
  framerate   = 10 ; FPS
  
  ; Set the plot dimensions and the peak colorabar scaling 
  ; Assume brightest pixel the first image in the sequence will set the top of the color bars
    Na_files = FILE_search('C:\IDL\Generic Model V2\read_write\'+write_directory+'\Na_Meteor_frame_*.fit')
    test_frame_Na = mrdfits(Na_files[0], 0, header, /silent) * brightness_multiplier_Na / 1.e3 
    dimensions    = size(test_frame_Na, /dimensions)
    xs            = dimensions[0]*4
    ys            = dimensions[1]*3

  ; Format layout and the time axis
    pos        = cgLayout([2,2], xgap = 0, ixmargin = 0, iymargin = 1, oymargin = 5, oxmargin = 0) 
    date_label = LABEL_DATE(DATE_FORMAT = ['%H:%I'])
    video      = IDLffVideoWrite(mpgFilename, format='mp4')
    stream     = video.AddVideoStream(xs, ys, framerate)
    !P.font    = 1  
    Simulated_UVVS_brightness_Na = MAKE_ARRAY(N_elements(DDR_Na), /float, VALUE = !values.f_nan) 
    loadgood

  for i = 0, N_elements(DDR_Na)-1 do begin
  
    ; only plot "key frames"
      junk = where(key_frames_Na eq i, count, /NULL)
      if count eq 0 then continue
    
    ; And only plate files that exist since some measurements are pre-impact  
      if not (FILE_TEST('C:\IDL\Generic Model V2\read_write\'+write_directory+'\Na_Meteor_frame_'+strcompress(string(i), /remove_all)+'.fit')) then continue

    ; We are going to create high-resolution PNG images, which we will add to the video stream.
    ; We create the high-resolution PNG image from PostScript intermediate files.
      cgPS_Open, Directory + Write_directory + '\movie.eps', /ENCAPSULATED, /Quiet
      cgDisplay, xs = xs, ys = ys
      device, SET_FONT = 'Helvetica Bold', /TT_FONT

      ; read the model images
        Meteor_frame_Na = mrdfits(Directory + write_directory + '\Na_Meteor_frame_'+strcompress(string(i), /remove_all)+'.fit', 0, Na_header, /silent) 
        platescale      = float(sxpar(Na_header, 'FOV')) / float(sxpar(Na_header, 'NAXIS1'))  ; "/pix image platescale
        
        if (N_elements(Meteor_frame_Na) eq 1) then continue
      
      ; read the UVVS pointing info  
        Pointing_Na   = mrdfits(Directory + write_directory + '\Na_Meteor_frame_'+strcompress(string(i), /remove_all)+'.fit', 2, /silent)
        frame_Na      =  brightness_multiplier_Na * Meteor_frame_Na / 1.e3     ; Convert brightness from R to KR

      ; Log the simulated brightness, assumed to be the mean of the brightness at the apertures 4 corner positions 
        Simulated_UVVS_brightness_Na[i]          = mean(interpolate(frame_Na, Pointing_Na.APERTURE_CORNERS[0,*], Pointing_Na.APERTURE_CORNERS[1,*], missing = !values.F_Nan), /NaN) $ ; transient
                                                   + Na_BG[i]                                                                                                                         ; background                                                                    ; background
        
        Range       = [-fov / (3600.*2.),  fov / (3600.*2.)]
        axis_format = { xticks:6, yticks:6, XRange:Range, xtitle:'Degrees', YRange:Range, charsize:0.75}
        number_of_colors = 8                                                   ; number of colors / labels on the color bar, loadgood has 8 
  
        P = pos[*,0]+[0.25, -.03, 0.25, .08] 
        cgImage, bytscl(alog10(frame_Na), Na_scale[0], Na_scale[1]), Position = P, AXKEYWORDS=axis_format, /Axes, /KEEP_ASPECT
        cgColorbar, Range = Na_scale, /Vertical, Position = [P[0]+.423, P[1], P[0]+.443, P[3]], $
          Title = Strcompress('Simulated Na Brightness Log!D10!N[kR]'), TLocation='Left', /right, divisions = number_of_colors, charsize = 0.75
        
        ; Mark the UVVS aperture projection, the coordinates here need to be in degrees  
          cgPolygon, [transpose(Pointing_Na.APERTURE_CORNERS[0,*]), Pointing_Na.APERTURE_CORNERS[0,0]]*platescale/3600.+range[0], $
                     [transpose(Pointing_Na.APERTURE_CORNERS[1,*]), Pointing_Na.APERTURE_CORNERS[1,0]]*platescale/3600.+range[0], COLOR='Snow', /fill
          
;          if keyword_set(debug) then begin
;            print_struct, DDR_Na_subset[i], ['TARGET_LOCAL_TIME', 'TARGET_ALTITUDE']
;            print, Pointing_Na.APERTURE_CORNERS*platescale/3600.+range[0]
;            print, 'UTC: ', Na_UTC_string[i]
;          endif

         ;----------------------------Generate the Latitude & Longitude grid to overlay on the column density----------------------------------
         ; Find the angular size
           CSPICE_SPKEZR, body, PDS_ET_Na[i], 'J2000', 'LT', viewpoint, BODY_state, ltime
           R_M = 206264.806 * atan(radii[0] / norm(BODY_state[0:2]))          ; Radius of the body in arcsec

         ; Find the sub-observer planetographic longitude and latitude
           cspice_subpnt, 'Near point: ellipsoid', Body, PDS_ET_Na[i], 'IAU_'+body, 'LT+S', viewpoint, $
             spoint, trgepc, srfvec
           f = ( radii[0]-radii[2] ) / radii[0]                               ; flatness parameter, mercury is round so this is zero
           cspice_recpgr, body, spoint, radii[0], f, spglon, spglat, spgalt
           SOP_lat = spglat * cspice_dpr()                                    ; planetographic latitude of the sub-observer point on the surface. (0 at equator, 90 at north pole, -90 at south pole)
           SOP_lon = spglon * cspice_dpr()                                    ; planetographic longitude of the sub-observer point on the surface. (0 to 360 west longitude).
           ;if keyword_set(debug) then print, 'Sub-observer Lat & Lon  =', SOP_lat, SOP_lon, ' Angular diameter of ', body, ' =', 2.*R_M, ' arcsec'

         ; setup the x-y grid image
           xdim     = float(sxpar(Na_header, 'NAXIS1'))
           ydim     = float(sxpar(Na_header, 'NAXIS2'))
           ctr_xpix = xdim / 2.
           ctr_ypix = ydim / 2.
           pix2km   = platescale * (1./R_M) * radii[0]                        ; "/pix * radii/" * km/radii
           x        = (dindgen(xdim) - ctr_xpix) * pix2km
           y        = (dindgen(xdim) - ctr_ypix) * pix2km
           xsq      = dblarr(xdim, ydim)                                      ; xsq is an image the size of the calib-img where the values are the horizontal distance from the "center" in km
           ysq      = dblarr(xdim, ydim)                                      ; ysq is an image the size of the calib-img where the values are the vertical distance from the center.
           for j = 0, xdim-1 do begin
             xsq[*,j] = x
             ysq[j,*] = y
           endfor

         ; Get the planetographic lon/lat at each pixel.
         ; This grid looks best if make the arrays big, rotate things to our "Sunward = -X, Observer-Planet = +Z" geometry with interpolation, then shrink them again.
           ysq = rebin(ysq, 1024, 1024)
           xsq = rebin(xsq, 1024, 1024)
           ob     = deprob(xsq, ysq, radii[0], radii[2], SOP_lat, SOP_lon)
           ob.lon[where(ob.lon eq -666, /Null)] = !values.F_nan
           ob.lat[where(ob.lat eq -666, /Null)] = !values.F_nan
           ob.lon = rot(ob.lon, sub_solar_lon*!radeg - SOP_lon - 90., /interp)
           ob.lat = rot(ob.lat, sub_solar_lon*!radeg - SOP_lon - 90., /interp)
           ob.lon = rebin(ob.lon, 1024, 1024)
           ob.lat = rebin(ob.lat, 1024, 1024)

         ; Mark the impact site's lat and lon
           cgcontour, ob.lat, /onimage, levels = float(STRMID(Surface_distribution, STRPOS(Surface_distribution, ',')+1, STRPOS(Surface_distribution, ']')-STRPOS(Surface_distribution, ',')-1)), $
             color = 'red', THICK = 1.5, LABEL = 0
           cgcontour, ob.lon, /onimage, levels = float(STRMID(Surface_distribution, STRPOS(Surface_distribution, '[')+1, STRPOS(Surface_distribution, ',')-STRPOS(Surface_distribution, '[')-1)), $
             color = 'red', THICK = 1.5, LABEL = 0, OLEVELS = OLEVELS
           ;cgcontour, ob.lon, /onimage, levels = sub_solar_lon*!radeg, color = 'YELLOW', THICK = 0.5 ; Mark the Sub-Solar Point

         ; Mark the lat & lon grid
           lat_contours = indgen(17)*10 - 80                                  ; every 10 deg lat
           lon_contours = indgen(12)*30                                       ; every 30 deg lon
           cgcontour, ob.lat, /onimage, levels = lat_contours, LABEL = 0, color = 'Snow', THICK = .5, CHARSIZE = .5
           cgcontour, ob.lon, /onimage, levels = lon_contours, LABEL = 1, color = 'Snow', THICK = .5, CHARSIZE = .5
           
           cgtext, -28, 25, 'Na', color = Na_Model_Color, charsize = 1.2
           
       ; Plot the time-series panel
         P = [pos[0,2]+0.065, pos[3,2]-.3, pos[2,3]-0.065, pos[3,2]+.03]
         cgplot, PDS_JD_Na, DDR_Na, /nodata, /xstyle, XTICKFORMAT='LABEL_DATE', XTICKUNITS = ['Minutes'], Xminor=5, XTICKINTERVAL=15, xtitle = strmid(Na_UTC_string[0], 0, 10)+' UTC', $
            Position=P, yr = [-1., 27.5], charsize = 1, /noerase, ERR_YHIGH = DDR_Na_err, ERR_Ylow = DDR_Na_err, /ERR_CLIP, ERR_WIDTH = 0.
       
        ; plot the rolling time ticker
           cgplot, [PDS_JD_Na[i], PDS_JD_Na[i]], [-10, 50], thick = 0.5, /overplot
          
        ; now the UVVS data
           cgplot, PDS_JD_Na, DDR_Na, color = Na_Data_Color, psym=16, XTICKUNITS = ['Minutes'], Xminor=5, XTICKINTERVAL=15, /xstyle, symsize = 0.5 , $
           Position=P, yr = [-1., 27.5], charsize = 1, /noerase, ERR_YHIGH = DDR_Na_err, ERR_Ylow = DDR_Na_err, /ERR_CLIP, ERR_WIDTH = 0., xtickformat = '(A1)

        ; plot the transient cloud
            plot_indices = where(finite(Simulated_UVVS_brightness_Na[0:i]), count_good_ind, /NULL)
            if count_good_ind gt 0 then cgplot, PDS_JD_Na[plot_indices], Simulated_UVVS_brightness_Na[plot_indices], color = 'Black', /overplot, thick = 4   
            ;cgplot, PDS_JD_Na[0:i], Simulated_UVVS_brightness_Na[0:i], color = 'Black', /overplot, thick = 4 

        ; plot the nominal background  
            cgplot, PDS_JD_Na[0:i], Na_BG[0:i], color = 'Black', /overplot, thick = 4, linestyle = 1                

        ; Annotate what's going on here  
          cgtext, PDS_JD_Na[20], 24.5, 'UVVS Na', color = Na_Data_Color, charsize = 1
          cgtext, mean(!x.crange), 24, 'Simulated Meteor Impact: ' + strmid(Meteor_impact_UTC, 11, 5) + ' at ' + $
            strmid(Surface_distribution, 7, 3) + cgsymbol('deg') + ' W Lon, '+ strmid(Surface_distribution, 11, 3) + cgsymbol('deg') + ' Lat', color = 'black', charsize = 1, alignment = .5
          cgtext, mean(!x.crange), 21, string(brightness_multiplier_Na * sxpar(Na_header, 'UPWARD_F') * 22.989769*1.e-3 / 6.02214e23, format = '(F4.2)') + ' kg Na, ' + Plume_Temperature + ' MBF', $
            color = 'black', charsize = 1, alignment = .5
          cgtext, PDS_JD_Na[0]-.01, 15, 'Brightness [kR]', orientation = 90, charsize = 1.

      cgPS_Close, width = xs, /PNG;, /Delete_PS ; Convert to PNG file.
      image = Read_PNG(Directory + write_directory + '\movie.png')
      image = image[*,0:xs-1,0:ys-1] ;not sure why this is needed but it is
      void = video -> Put(stream, image)

  endfor
  video -> Cleanup

  ; Compare the results by correlation
    keep_Na   = where(finite(Simulated_UVVS_brightness_Na))
    correl_Na = correlate(Simulated_UVVS_brightness_Na[keep_Na], DDR_Na[keep_Na])
    Na_Chisq  = total( ( (DDR_Na - Simulated_UVVS_brightness_Na) / DDR_Na_err)^2., /NAN) / total(finite(Simulated_UVVS_brightness_Na))
    
    Print, 'Correlation of results Na[Key_frames]:', correl_Na
    Print, 'Chi Squared of results Na[Key_frames]:', Na_Chisq
    save, correl_Na, filename = Directory+ 'Correlation_results\' + write_directory + '_correlation_Event2_Round_0.sav'
    writecol, Directory + write_directory + '\Na_transient.txt', Na_UTC_string, Simulated_UVVS_brightness_Na
    ;stop
end