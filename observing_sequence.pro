Pro Observing_Sequence, Meteor_impact_UTC = Meteor_impact_UTC, Plume_Temperature = Plume_Temperature, $
                        Surface_distribution = Surface_distribution, loop_times = loop_times

; Runtimes:

; Background: 1.e26 1200K MBF, 100 loops, 0.166 Days duration = 2 days runtime
; Meteor:     5.e27 3100K MBF, 12 loops = 12.5 hours runtime

; Most Promising Results so far: 
;      Meteor_impact_UTC             = '2011-08-04 01:40:00'      ; time of the impact
;      Plume_Temperature             = '2000K'                    ; temperature of the impact vapour
;      Surface_distribution = 'Point_[0, -30]'                    ; Location of the impactor
;      loop_times                    = 40.                        ; Bear minimum for any reasonable S/N
;      Na_Lofted                     = 7.e28                      ; seems like a lot

; Round 2
      Meteor_impact_UTC             = '2011-08-04 01:40:00'      ; time of the impact
      Plume_Temperature             = '3500K'                    ; temperature of the impact vapour
      Surface_distribution = 'Point_[345, -45]'                    ; Location of the impactor
      loop_times                     = 20.                       ; Bear minimum for any reasonable S/N
      Na_Lofted                     = 1.4e28                      ; seems like a lot
      Mg_Lofted                     = 5.6e28                      ; seems like a lot
      Brightness_multiplier_Na      = 1.
      Brightness_multiplier_Mg      = 1.

COMMON Output_shared, Plot_range, Output_Size_In_Pixels, Output_Title, Center_in_frame, viewpoint, FOV, N_ticks, Tickstep, Observatory, Above_Ecliptic, Boresight_Pixel
COMMON Model_shared, Body, Ephemeris_time, Seed, Directory, Particle_data, Line_data, Debug

;Brightness_multiplier_Na      = 11.2
;Brightness_multiplier_Mg      = 24.  
;Na_Lofted                     = 1.e27
;Mg_Lofted                     = 4.*Na_Lofted

;key_frames                    = [7, 15, 28, 48]                ; match these 4 and you've done it!
;key_frames                    = [5, 7, 16, 21, 27, 45, 48, 51]  ; match these 8 and you've done it!
key_frames                    = indgen(29) * 2 
;key_frames                    = indgen(28) * 2 + 2

; Calculate the flatness coefficient and planetary radius
  cspice_bodvrd, body, 'RADII', 3, radii
  re   = radii[0]
  rp   = radii[2]
  flat = ( re - rp ) / re

; Print the longitude and latitude of equatorial dawn at the impact time
  cspice_utc2et, Meteor_impact_UTC, meteor_ET
  cspice_subslr, 'Near point: ellipsoid', body, meteor_ET, 'IAU_'+body, 'none', viewpoint, sub_solar_point_planet_frame, trgepc, srfvec ; not sure why we'd need light time
  cspice_recpgr, body, sub_solar_point_planet_frame, re, flat, sub_solar_lon, sub_solar_lat, sub_solar_radius ; Convert to planetographic latitude and longitude, in the 'IAU_' body fixed frame  
  print, 'Coordinates of Equatorial Dawn [Lon, Lat]', (sub_solar_lon + !pi/2. + 2.*!pi) mod (2.*!pi)*!radeg, sub_solar_lat*!radeg

; Load UVVS DDR
; This is downloaded from: https://pds-geosciences.wustl.edu/messenger/mess-e_v_h-mascs-3-virs-cdr-caldata-v1/messmas_2101/data/ddr/ob2/uvvs_atmosphere/07/
  UVVS_DDR_Na      = read_mascs_ddr(Directory+'MESSENGER_UVVS\ud_02_ns_na.dat') 
  UVVS_UTC_TIME_Na = string(UVVS_DDR_Na.UTC_TIME)  ; convert time from byte to string YYDOYTHH:MM:SS.00
  UVVS_DDR_Mg      = read_mascs_ddr(Directory+'MESSENGER_UVVS\ud_02_ns_mg.dat') 
  UVVS_UTC_TIME_Mg = string(UVVS_DDR_Mg.UTC_TIME)  ; convert time from byte to string YYDOYTHH:MM:SS.00

; Load Tim's save file
  restore, Directory+'MESSENGER_UVVS\orbit 278 radiance for Carl.sav', /verbose
  Tim_utc_time = utc_time
  Tim_RADIANCE_KR = RADIANCE_KR

; Check for a match between Tim's save file and the PDS
  start  = where(UVVS_UTC_TIME_Na eq Tim_utc_time[0])
  stop   = where(UVVS_UTC_TIME_Na eq Tim_utc_time[-1])
  
  date   = Tim_utc_time
  Tim_ET = dblarr(N_elements(Tim_utc_time))
  JD     = dblarr(N_elements(Tim_utc_time))
  year = '20' + StrMid(StrTrim(date,2), 0, 2)
  dayofyear = StrMid(StrTrim(date,2), 2, 3)
  time = StrMid(StrTrim(date,2), 6, 11)
  CALDAT, JULDAY(1, dayofyear, year), month, day
  UTC_string = year + '-' + string(month ,format = "(I2.2)") + '-' + string(day ,format = "(I2.2)") + ' ' + time
  for i = 0, N_elements(Tim_utc_time)-1 do begin
    cspice_utc2et, UTC_string[i], ET
    Tim_ET[i] = ET
    cspice_et2utc, et, 'J', 6, JD_string
    JD[i] = strmid(JD_string, 3)
  endfor
  
  date = UVVS_UTC_TIME_Na[Start:Stop]
  PDS_ET = dblarr(N_elements(date))
  PDS_JD_Na = dblarr(N_elements(date))
  year = '20' + StrMid(StrTrim(date,2), 0, 2)
  dayofyear = StrMid(StrTrim(date,2), 2, 3)
  time = StrMid(StrTrim(date,2), 6, 11)
  CALDAT, JULDAY(1, dayofyear, year), month, day
  PDS_UTC_string = year + '-' + string(month ,format = "(I2.2)") + '-' + string(day ,format = "(I2.2)") + ' ' + time
  for i = 0, N_elements(date)-1 do begin
    cspice_utc2et, PDS_UTC_string[i], ET
    PDS_ET[i] = ET
    cspice_et2utc, et, 'J', 6, JD_string
    PDS_JD_Na[i] = strmid(JD_string, 3)
  endfor

  date = UVVS_UTC_TIME_Mg[1989:2396]
  PDS_ET = dblarr(N_elements(date))
  PDS_JD_Mg = dblarr(N_elements(date))
  year = '20' + StrMid(StrTrim(date,2), 0, 2)
  dayofyear = StrMid(StrTrim(date,2), 2, 3)
  time = StrMid(StrTrim(date,2), 6, 11)
  CALDAT, JULDAY(1, dayofyear, year), month, day
  PDS_UTC_string = year + '-' + string(month ,format = "(I2.2)") + '-' + string(day ,format = "(I2.2)") + ' ' + time
  for i = 0, N_elements(date)-1 do begin
    cspice_utc2et, PDS_UTC_string[i], ET
    PDS_ET[i] = ET
    cspice_et2utc, et, 'J', 6, JD_string
    PDS_JD_Mg[i] = strmid(JD_string, 3)
  endfor

  ; get a subset of the PDS data with a SNR threshold
    SNR_threshold = 3.
    DDR_Na    = UVVS_DDR_Na[Start:Stop].TOTAL_RADIANCE_KR
    keep_Na   = where(UVVS_DDR_Na[Start:Stop].TOTAL_RADIANCE_SNR gt SNR_threshold, /Null)
    DDR_Na    = DDR_Na[keep_Na]
    PDS_JD_Na = PDS_JD_Na[Keep_Na]
    
    DDR_Mg    = UVVS_DDR_Mg[1989:2396].TOTAL_RADIANCE_KR
    keep_Mg   = where(UVVS_DDR_Mg[1989:2396].TOTAL_RADIANCE_SNR gt SNR_threshold, /Null)
    DDR_Mg    = DDR_Mg[keep_Mg] 
    PDS_JD_Mg = PDS_JD_Mg[Keep_Mg]
  
;window, 0, xs = 1800, ys = 800
;  
;    cgplot, PDS_JD_Na, DDR_Na, color = 'grey', psym=16, /xstyle, ytitle = 'Na D Brightness [KR]', XTICKUNITS = ['Time'], $
;    XTICKFORMAT='LABEL_DATE', xtitle = strmid(utc_string[0], 0, 10)+' UTC', Position=P, yr = [-1., 18.], charsize = 1, /noerase
;    cgplot, PDS_JD_Mg, DDR_Mg, color = 'red', psym=16, /overplot

days_since_meteor_impact = (Tim_ET - meteor_ET) / (24.*3600.)

; Format these inputs into a string and make a new write direcctory for this observing sequence  
  cspice_et2utc, meteor_ET, 'ISOC', 0, impact_UTC
  Surface_distribution_string = strmid(Surface_distribution, 7)
  ;Surface_distribution_string = repstr(Surface_distribution_string,', ','lon')
  ;Surface_distribution_string = repstr(Surface_distribution_string)
  Surface_distribution_string = strcompress(repstr(Surface_distribution_string,']','lat'), /remove_all)
  impact_UTC                  = repstr(impact_UTC,':','-') ; need to replace the colons in the directory string, use '-'
  
  cd, directory
  write_directory = impact_UTC + '---' + Surface_distribution_string + '---' + Plume_Temperature
  FILE_MKDIR, write_directory, /NOEXPAND_PATH

  ;Profiler, /RESET
  ;profiler
  ;profiler, /system

restore, Directory+'boresight_pixels.sav'
  boresight_pixels = intarr(2, N_elements(Tim_utc_time)) ;log the boresight pixel location at each time
  for i = 0, N_elements(Tim_utc_time)-1 do begin
    junk = where(key_frames eq i, count, /NULL)
    if count eq 0 then continue
    if days_since_meteor_impact[i] lt 0. then continue ; Can't model an event that hasn't happened yet    
    if FILE_TEST('C:\IDL\Generic Model V2\read_write\'+write_directory+'\Na_Meteor_frame_'+strcompress(string(i), /remove_all)+'.fit') then begin ; skip if we've done this before
      restore, Directory+'boresight_pixels.sav'
      continue 
    endif
;    
;     Generic_Model, Time_range_this_run = days_since_meteor_impact[i], Output_title_this_run = write_directory + '\Na_Meteor_frame_'+strcompress(string(i), /remove_all), $
;                    test_particle_this_run = 'Na', Line_this_run = 'Na-D', UTC_this_run = UTC_string[i], $
;                    Surface_distribution_this_run = Surface_distribution, Speed_distribution_this_run = 'MBF_'+plume_temperature, $
;                    loop_times_this_run = loop_times, Upward_flux_at_exobase_this_run = Na_Lofted                
;     Generic_Model, Time_range_this_run = days_since_meteor_impact[i], Output_title_this_run = write_directory + '\Mg_Meteor_frame_'+strcompress(string(i), /remove_all), $
;                    test_particle_this_run = 'Mg', Line_this_run = 'Mg-2853',  UTC_this_run = UTC_string[i], $
;                    Surface_distribution_this_run = Surface_distribution, Speed_distribution_this_run = 'MBF_'+plume_temperature, $
;                    loop_times_this_run = loop_times, Upward_flux_at_exobase_this_run = Mg_Lofted
;                    
      Generic_Model, Time_range_this_run = days_since_meteor_impact[i], Output_title_this_run = write_directory + '\Na_Meteor_frame_'+strcompress(string(i), /remove_all), $
        test_particle_this_run = 'Na', Line_this_run = 'Na-D', UTC_this_run = UTC_string[i], $
        Surface_distribution_this_run = Surface_distribution, Speed_distribution_this_run = 'Cone_'+plume_temperature, $
        loop_times_this_run = loop_times, Upward_flux_at_exobase_this_run = Na_Lofted
      
      Generic_Model, Time_range_this_run = days_since_meteor_impact[i], Output_title_this_run = write_directory + '\Mg_Meteor_frame_'+strcompress(string(i), /remove_all), $
        test_particle_this_run = 'Mg', Line_this_run = 'Mg-2853',  UTC_this_run = UTC_string[i], $
        Surface_distribution_this_run = Surface_distribution, Speed_distribution_this_run = 'Cone_'+plume_temperature, $
        loop_times_this_run = loop_times, Upward_flux_at_exobase_this_run = Mg_Lofted
                    
    boresight_pixels[*,i] = boresight_pixel
  endfor

;  profiler, /report, data = data
;  profiler, /clear
;  sorted = sort(-data.time)
;  print, data[sorted[0:9]], format = '(A-20, I7, F12.5, F10.5, I9)'

  ;-----------------------MOVIE-----------------------------------------
  ;files       = file_search(Directory+write_directory+'\Na*frame_*.fit', count = N_frames) ; this is just used to count frames, since the returned array has files listed out of order
  mpgFilename = Directory+Write_directory+'\test_cone.mp4'
  framerate   = 10 ; FPS
  
  ; get the size and set the peak ccolorabar scaling 
  ; (assume brightest pixel the first image in the sequence will set the top of the color bars)
    test_frame_Na  = mrdfits(Directory+write_directory+'\Na_Meteor_frame_'+strcompress(string(key_frames[0]), /remove_all)+'.fit', 0, header, /silent) * brightness_multiplier_Na / 1.e3 
    test_frame_mg  = mrdfits(Directory+write_directory+'\Mg_Meteor_frame_'+strcompress(string(key_frames[0]), /remove_all)+'.fit', 0, header, /silent) * brightness_multiplier_Mg 
    dimensions  = size(test_frame_Na, /dimensions)
    xs          = dimensions[0]*8
    ys          = dimensions[1]*6
 
  ; Format layout and the time axis
    pos        = cgLayout([2,2], xgap = 0, ixmargin = 0, iymargin = 1, oymargin = 5, oxmargin = 0) 
    date_label = LABEL_DATE(DATE_FORMAT = ['%H:%I'])
    video      = IDLffVideoWrite(mpgFilename, format='mp4')
    stream     = video.AddVideoStream(xs, ys, framerate)
    !P.font    = 1  
;    MinValue   = 0
;    MaxValue   = 500.
;    MinValue_CD= 1.e0
;    MaxValue_CD= 5000.  
    Simulated_UVVS_brightness_Na = MAKE_ARRAY(N_elements(Tim_utc_time), /float, VALUE = !values.f_nan) 
    Simulated_UVVS_brightness_Mg = MAKE_ARRAY(N_elements(Tim_utc_time), /float, VALUE = !values.f_nan) 
    loadgood

  for i  = 0, N_elements(Tim_utc_time) -1 do begin
    
    junk = where(key_frames eq i, count, /NULL)
    if count eq 0 then continue

    ; We are going to create high-resolution PNG images, which we will add
    ; to the video stream. We create the high-resolution PNG image from PostScript
    ; intermediate files.
      cgPS_Open, Directory + Write_directory + '\movie.eps', /ENCAPSULATED, /Quiet
      cgDisplay, xs = xs, ys = ys
      device, SET_FONT = 'Helvetica Bold', /TT_FONT

        Meteor_frame_Na        = mrdfits(Directory + write_directory + '\Na_Meteor_frame_'+strcompress(string(i), /remove_all)+'.fit', 0, Na_header, /silent)
        Meteor_frame_Mg        = mrdfits(Directory + write_directory + '\Mg_Meteor_frame_'+strcompress(string(i), /remove_all)+'.fit', 0, Mg_Header, /silent)
        
        if (N_elements(Meteor_frame_Na) eq 1) or (N_elements(Meteor_frame_Mg) eq 1) then continue

        ;frame       = (background_frame + Meteor_frame) / 1.e3 ; Convert brightness from R to KR
        frame_Na      =  brightness_multiplier_Na * Meteor_frame_Na / 1.e3 ; Convert brightness from R to KR
        frame_Mg      =  brightness_multiplier_Mg * Meteor_frame_Mg 
      
      ; Log the brightness and Mark the boresight's tanget point  
        Simulated_UVVS_brightness_Na[i] = mean( frame_Na[ boresight_pixels[0,i]-1:boresight_pixels[0,i]+1, boresight_pixels[1,i]-1:boresight_pixels[1,i]+1 ] ) 
        Simulated_UVVS_brightness_Mg[i] = mean( frame_Mg[ boresight_pixels[0,i]-1:boresight_pixels[0,i]+1, boresight_pixels[1,i]-1:boresight_pixels[1,i]+1 ] ) 
        frame_Na[ boresight_pixels[0,i]-1:boresight_pixels[0,i]+1, boresight_pixels[1,i]-1:boresight_pixels[1,i]+1 ] = max(test_frame_Na)    
        frame_Mg[ boresight_pixels[0,i]-1:boresight_pixels[0,i]+1, boresight_pixels[1,i]-1:boresight_pixels[1,i]+1 ] = max(test_frame_Mg)
  
        Range       = [-fov / (3600.*2.),  fov / (3600.*2.)]
        axis_format = { xticks:6, yticks:6, XRange:Range, xtitle:'Degrees', YRange:Range, charsize:0.75}
  
        number_of_colors = 8                                                     ;number of colors / labels on the color bar, loadgood has 8 
  
        P = pos[*,0]+[0, -.03, 0, .08] 
        cgImage, bytscl(alog10(frame_Na), 0., max(alog10(test_frame_Na))), Position = P, AXKEYWORDS=axis_format, /Axes, /KEEP_ASPECT
        cgColorbar, Range=[0., max(alog10(test_frame_Na))], /Vertical, Position = [.445, P[1], .465, P[3]], $
          Title = Strcompress('Simulated Na Brightness Log!D10!N[KR]'), TLocation='Left', /right, divisions = number_of_colors, charsize = 0.75

        P = pos[*,1]+[0.015, -.03, 0., .08]
        cgImage, bytscl(alog10(frame_Mg), 2., max(alog10(test_frame_Mg))), AXKEYWORDS=axis_format, /Axes, /KEEP_ASPECT, Position = P, /noerase
        cgColorbar, Range=[2., max(alog10(test_frame_Mg), /Nan)], /Vertical, Position = [.93, P[1], .95, P[3]], $
          Title = Strcompress('Simulated Mg Brightness Log!D10!N[R]'), TLocation='Left', /right, divisions = number_of_colors, charsize = 0.75, $
          /noerase
          
            ;----------------------------Generate the Latitude & Longitude grid to overlay on the column density----------------------------------
            ; Find the position angle of the north pole vector, that is, the separation between Body's north pole and celestial north pole.
            CSPICE_SPKEZR, body, Tim_ET[i], 'J2000', 'LT', viewpoint, BODY_state, ltime
            CSPICE_RECRAD, BODY_state[0:2], dist, ra, dec ;Convert rectangular coordinates to RA and Dec
            R_M = 206264.806 * atan(radii[0] / norm(BODY_state[0:2]))          ; Radius of the body in arcsec
    
            ; Find a vector from body's center to it's north pole in the J2000 frame at the ephemeris time, NEGLECT the planet's oblateness
            North_body_fixed = [0.,0.,1.]*radii[2]                             ; location of the North Pole in body-fixed coords
            cspice_pxform, 'IAU_'+Body, 'J2000', Tim_ET[i] - ltime, To_J2000   ; Find body-fixed coords to J2000 rotation matrix
            cspice_mxv, To_J2000, North_body_fixed, North_J2000                ; Rotate to J2000
    
            ; Get the RA and Dec of the the body's north pole
            Pole_state = body_state + North_J2000
            CSPICE_RECRAD, Pole_state, dist, Pole_ra, Pole_dec ;Convert rectangular coordinates to RA and Dec
    
            ; Precess each RA and Dec to the Current Epoch. First we'll find the year and fraction of a year. we're interested in:
            Current_epoch = 2000. + Tim_ET[i] / 3.1556926e7                    ; J2000 + seconds past J2000 / seconds per year
            precess, Pole_ra, Pole_DEC, 2000, Current_epoch, /RADIAN           ; Precessed Ra and dec is now applied to Pole location
            precess, ra, DEC, 2000, Current_epoch, /RADIAN                     ; Precessed Ra and dec is now applied to body center location
            Delta_RA      = Pole_ra - ra
            Delta_Dec     = Pole_Dec - dec
            theta         = sqrt((Delta_RA*cos(Dec))^2.D + Delta_Dec^2.D)
            NPAng         = signum(Delta_RA)*Acos(Delta_Dec / theta) / cspice_rpd()
    
            ; Find the sub-observer planetographic longitude and latitude
            cspice_subpnt, 'Near point: ellipsoid', Body, Tim_ET[i], 'IAU_'+body, 'LT+S', viewpoint, $
              spoint, trgepc, srfvec
            f = ( radii[0]-radii[2] ) / radii[0]                               ; flatness parameter, mercury is round so this is zero
            cspice_recpgr, body, spoint, radii[0], f, spglon, spglat, spgalt
            lat_se = spglat * cspice_dpr()                                     ; planetographic latitude of the sub-observer point on the surface. (0 at equator, 90 at north pole, -90 at south pole)
            lon_se = spglon * cspice_dpr()                                     ; planetographic longitude of the sub-observer point on the surface. (0 to 360 west longitude).
            if keyword_set(debug) then print, 'North Pole Position Angle (CCW, E of N) = ', NPAng, $
              ' sub-observer Lat & Lon  =', lat_se, lon_se, ' Angular diameter of ', body, ' =', 2.*R_M, ' arcsec'
    
            ; setup the x-y grid image
            platescale = 90.*3600./128.                                        ; "/pix HACK: hard code this using the headers
            xdim     = 128.
            ydim     = 128.
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
    
            ; get the planetographic lon/lat at each pixel
            ob           = deprob(xsq, ysq, radii[0], radii[2], lat_se, lon_se, npang=NPAng)
            ob.lon[where(ob.lon eq -666, /Null)] = !values.F_nan
            ob.lat[where(ob.lat eq -666, /Null)] = !values.F_nan
            ;ob.lon[where(frame gt 0., /Null)] = !values.F_nan
            ;ob.lat[where(frame gt 0., /Null)] = !values.F_nan
    
            lat_contours = indgen(17)*10 - 80                                  ; every 10 deg lat
            lon_contours = indgen(24)*15                                       ; every 15 deg lon
            
            cgcontour, ob.lat, /onimage, levels = lat_contours, LABEL = 0, color = 'green', THICK = .5;, C_Spacing = 5
            ;cgcontour, ob.lon, /onimage, levels = lon_contours, LABEL = 0, color = 'green', THICK = .5;, C_Spacing = 5 ;'snow'
           
             cgcontour, ob.lon, /onimage, NLEVELS = 10, LABEL = 0, color = 'green', THICK = .5
            
            cgcontour, ob.lat, /onimage, levels = float(STRMID(Surface_distribution, STRPOS(Surface_distribution, ',')+1, STRPOS(Surface_distribution, ']')-STRPOS(Surface_distribution, ',')-1)), $
              LABEL = 0, color = 'red', THICK = .5
            ;cgcontour, ob.lon, /onimage, levels = float(STRMID(Surface_distribution, STRPOS(Surface_distribution, '[')+1, STRPOS(Surface_distribution, ',')-STRPOS(Surface_distribution, '[')-1)), $
            ;  LABEL = 0, color = 'red', THICK = .5
            ;cgcontour, ob.lon, /onimage, levels = sub_solar_lon*!radeg, color = 'red'
          
        P = [pos[0,2]+0.065, pos[1,2], pos[2,3]-0.065, pos[3,2]]
        cgplot, PDS_JD_Na, DDR_Na, color = 'Black', psym=16, /xstyle, YSTYLE=9, ytitle = 'Na D Brightness [KR]', XTICKUNITS = ['Time'], $
          XTICKFORMAT='LABEL_DATE', xtitle = strmid(utc_string[0], 0, 10)+' UTC', Position=P, yr = [-1., 18.], charsize = 1, /noerase
        cgplot, PDS_JD_Mg, DDR_Mg*10., color = 'red', psym=16, XRANGE = !X.CRANGE, Position=P, charsize = 1, /noerase, ystyle = 16, XTICKUNITS = ['Time'], $
          XTICKFORMAT='LABEL_DATE', yr = [-1., 18.]
        cgAXIS, YAXIS=1, YRANGE = !Y.CRANGE*100., YSTYLE = 1, YTITLE = 'Mg 2853A Brightness [R]', charsize = 1., color = 'red'
           
        ;cgplot, JD, Tim_RADIANCE_KR, color = 'Black', psym=16, /overplot
        cgplot, JD[0:i], Simulated_UVVS_brightness_Na[0:i], color = 'grey', psym=14, /overplot, symsize = 0.75      
        cgplot, JD[0:i], Simulated_UVVS_brightness_Mg[0:i]/100., color = 'pink', psym=14, /overplot, symsize = 0.75   

        cgplot, [JD[i], JD[i]], [0, 20], /overplot
        cgtext, JD[87],15, 'UVVS Na (SNR > 3)', color = 'Black', charsize = 1
        cgtext, JD[87],13, 'UVVS Mg (SNR > 3)', Color = 'red', charsize = 1
        cgtext, JD[87],11, 'Simulation: ' + strmid(Meteor_impact_UTC, 10, 6) + ' [lon, lat] = ' + strmid(Surface_distribution, 6), color = 'black', charsize = 1
        cgtext, JD[87] + .72/24., 9, Plume_Temperature + ' MBF', color = 'black', charsize = 1
        cgtext, JD[87] + .70/24., 7, strcompress( round(brightness_multiplier_Na * sxpar(Na_header, 'UPWARD_F') * 22.989769*1.e-3 / 6.02214e23)) + ' kg (Na)', color = 'grey', charsize = 1
        cgtext, JD[87] + 1.5/24., 7, strcompress( round(brightness_multiplier_Mg * sxpar(Mg_header, 'UPWARD_F') * 24.305*1.e-3 / 6.02214e23)) + ' kg (Mg)', color = 'pink', charsize = 1

      cgPS_Close, width = xs, /PNG;, /Delete_PS ; Convert to PNG file.

      image = Read_PNG(Directory + write_directory + '\movie.png')
      image = image[*,0:xs-1,0:ys-1] ;not sure why this is needed but it is
      ;File_Delete, Directory + write_directory + '\movie.png'
      void = video -> Put(stream, image)
  endfor
  video -> Cleanup
  
  ; Compare the results by correlation
    correl = correlate(Simulated_UVVS_brightness_Na[key_frames], Tim_RADIANCE_KR[key_frames])
    Print, 'correlation of Na results [Key_frames]]:', correl
    save, correl, filename = Directory+ 'Correlation_results\' + write_directory + '_Na_correlation_Round_Hi_res_cone.sav'

end