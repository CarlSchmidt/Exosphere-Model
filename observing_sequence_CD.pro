Pro Observing_Sequence

; Runtimes:

; Background: 1.e26 1200K MBF, 100 loops, 0.166 Days duration = 2 days runtime
; Meteor:     5.e27 3100K MBF, 12 loops = 12.5 hours runtime

; Promising Results: 

;  Meteor_impact_UTC             = '2011-08-04 02:00:00'      ; time of the impact
;  Plume_Temperature             = '3000K'                    ; temperature of the impact vapour
;  Surface_distribution_this_run = 'Point_[50, -20.0]'        ; Location of the impactor
;  loop_times_this_run           = 30.
;  Upward_flux_at_exobase_this_run = 5.e27

COMMON Output_shared, Plot_range, Output_Size_In_Pixels, Output_Title, Center_in_frame, viewpoint, FOV, N_ticks, Tickstep, Observatory, Above_Ecliptic, Boresight_Pixel
COMMON Model_shared, Body, Ephemeris_time, Seed, Directory, Particle_data, Line_data, Debug

Meteor_impact_UTC             = '2011-08-04 02:15:00'      ; time of the impact
Plume_Temperature             = '6000K'                    ; temperature of the impact vapour
Surface_distribution_this_run = 'Point_[20, 20]'           ; Lon, Lat location of the impactor (Equatorial Dawn is 30,0)
Brightness_multiplier         = 10.  
loop_times_this_run           = 2. 
Na_Lofted                     = 1.e27
Mg_Lofted                     = 5.e27
key_frames                    = [6, 15, 28, 45]           ; match these 4 and you've done it!
;key_frames                    = [6, 15, 30, 45]            ; match these 4 and you've done it!

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

  date = UVVS_UTC_TIME_Mg[Start:Stop]
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
    DDR_Na  = UVVS_DDR_Na[Start:Stop].TOTAL_RADIANCE_KR
    keep_Na = where(UVVS_DDR_Na[Start:Stop].TOTAL_RADIANCE_SNR gt SNR_threshold, /Null)
    DDR_Na  = DDR_Na[keep_Na]
    DDR_Na  = PDS_JD_Na[Keep_Na]
    DDR_Mg    = UVVS_DDR_Mg[Start:Stop].TOTAL_RADIANCE_KR
    keep_Mg   = where(UVVS_DDR_Mg[Start:Stop].TOTAL_RADIANCE_SNR gt SNR_threshold, /Null)
    DDR_Mg    = DDR_Mg[keep_Mg] 
    PDS_JD_Mg = PDS_JD_Mg[Keep_Mg]
  
;window, 0, xs = 1800, ys = 800
;  
;    cgplot, PDS_JD_Na, DDR_Na, color = 'grey', psym=16, /xstyle, ytitle = 'Na D Brightness [KR]', XTICKUNITS = ['Time'], $
;    XTICKFORMAT='LABEL_DATE', xtitle = strmid(utc_string[0], 0, 10)+' UTC', Position=P, yr = [-1., 18.], charsize = 1, /noerase
;    cgplot, PDS_JD_Mg, DDR_Mg, color = 'red', psym=16, /overplot
;
;stop
;

days_since_meteor_impact = (Tim_ET - meteor_ET) / (24.*3600.)

; Format these inputs into a string and make a new write direcctory for this observing sequence  
  cspice_et2utc, meteor_ET, 'ISOC', 0, impact_UTC
  Surface_distribution_string = strmid(Surface_distribution_this_run, 7)
  Surface_distribution_string = repstr(Surface_distribution_string,', ','lon')
  Surface_distribution_string = repstr(Surface_distribution_string,']','lat')
  impact_UTC                  = repstr(impact_UTC,':','-') ; need to replace the colons in the directory string, use '-'
  cd, directory
  write_directory = impact_UTC + '---' + Surface_distribution_string + '---' + Plume_Temperature
  FILE_MKDIR, write_directory, /NOEXPAND_PATH

  ;Profiler, /RESET
  ;profiler
  ;profiler, /system

restore, Directory+'boresight_pixels.sav'
;  boresight_pixels = intarr(2, N_elements(Tim_utc_time)) ;log the boresight pixel location at each time
;  for i = 0, N_elements(Tim_utc_time)-1 do begin
;    junk = where(key_frames eq i, count, /NULL)
;    if count eq 0 then continue
;    ;if (i gt 50) then continue
;    if days_since_meteor_impact[i] lt 0. then continue ; Can't model an event that hasn't happened yet
;     ;Generic_Model, Output_title_this_run = 'Background_frame_'+strcompress(string(i), /remove_all), UTC_this_run = UTC_string[i]
;     
;     Generic_Model, Time_range_this_run = days_since_meteor_impact[i], Output_title_this_run = write_directory + '\Mg_Meteor_frame_'+strcompress(string(i), /remove_all), $
;                    test_particle_this_run = 'Mg', Line_this_run = 'Mg-2853',  UTC_this_run = UTC_string[i], $
;                    Surface_distribution_this_run = Surface_distribution_this_run, Speed_distribution_this_run = 'MBF_'+plume_temperature, $
;                    loop_times_this_run = loop_times_this_run, Upward_flux_at_exobase_this_run = Mg_Lofted
;     
;     Generic_Model, Time_range_this_run = days_since_meteor_impact[i], Output_title_this_run = write_directory + '\Na_Meteor_frame_'+strcompress(string(i), /remove_all), $
;                    test_particle_this_run = 'Na', Line_this_run = 'Na-D', UTC_this_run = UTC_string[i], $
;                    Surface_distribution_this_run = Surface_distribution_this_run, Speed_distribution_this_run = 'MBF_'+plume_temperature, $
;                    loop_times_this_run = loop_times_this_run, Upward_flux_at_exobase_this_run = Na_Lofted        
;                    
;    boresight_pixels[*,i] = boresight_pixel
;    ;if i eq 5 then break
;  endfor

;  profiler, /report, data = data
;  profiler, /clear
;  sorted = sort(-data.time)
;  print, data[sorted[0:9]], format = '(A-20, I7, F12.5, F10.5, I9)'

  ;-----------------------MOVIE-----------------------------------------
  files       = file_search(Directory+write_directory+'\*frame_*.fit', count = N_frames) ; this is just used to count frames, since the returned array has files listed out of order
  mpgFilename = Directory+Write_directory+'\test.mp4'
  framerate   = 10 ; FPS
  
  ; get the size and set the peak ccolorabar scaling 
  ; (assume brightest pixel the first image in the sequence will set the top of the color bars)
    test_frame     = mrdfits(files[0], 0, header, /silent) * brightness_multiplier / 1.e3 
    test_frame_CD  = mrdfits(files[0], 1, header, /silent) * brightness_multiplier
    dimensions  = size(test_frame, /dimensions)
    xs          = dimensions[0]*8
    ys          = dimensions[1]*6
 
  ; Format layout and the time axis
    pos        = cgLayout([2,2], xgap = 0, ixmargin = 0, iymargin = 0, oymargin = 4, oxmargin = 0) 
    date_label = LABEL_DATE(DATE_FORMAT = ['%H:%I'])
    video      = IDLffVideoWrite(mpgFilename, format='mp4')
    stream     = video.AddVideoStream(xs, ys, framerate)
    !P.font    = 1  
;    MinValue   = 0
;    MaxValue   = 500.
;    MinValue_CD= 1.e0
;    MaxValue_CD= 5000.  
    Simulated_UVVS_brightness = MAKE_ARRAY(N_elements(Tim_utc_time), /float, VALUE = !values.f_nan) ;Initialize 
    loadgood
 
  for i  = 0, N_elements(Tim_utc_time) -1 do begin
    
    junk = where(key_frames eq i, count, /NULL)
    if count eq 0 then continue
    ;if (i gt 50) then continue
    
    ; We are going to create high-resolution PNG images, which we will add
    ; to the video stream. We create the high-resolution PNG image from PostScript
    ; intermediate files.
      cgPS_Open, Directory + Write_directory + '\movie.eps', /ENCAPSULATED, /Quiet
      cgDisplay, xs = xs, ys = ys
      device, SET_FONT = 'Helvetica Bold', /TT_FONT

        Meteor_frame           = mrdfits(Directory + write_directory + '\Mg_Meteor_frame_'+strcompress(string(i), /remove_all)+'.fit', 0, header, /silent)
        Meteor_frame_CD        = mrdfits(Directory + write_directory + '\Mg_Meteor_frame_'+strcompress(string(i), /remove_all)+'.fit', 1, /silent)

        ;frame       = (background_frame + Meteor_frame) / 1.e3 ; Convert brightness from R to KR
        frame       =  brightness_multiplier * Meteor_frame / 1.e3 
        frame_CD    =  brightness_multiplier * Meteor_frame_CD 
      
      ; Log the brightness and Mark the boresight's tanget point  
        Simulated_UVVS_brightness[i] = mean( frame[ boresight_pixels[0,i]-1:boresight_pixels[0,i]+1, boresight_pixels[1,i]-1:boresight_pixels[1,i]+1 ] ) 
        frame[ boresight_pixels[0,i]-1:boresight_pixels[0,i]+1, boresight_pixels[1,i]-1:boresight_pixels[1,i]+1 ]    = max(test_frame)    
        frame_CD[ boresight_pixels[0,i]-1:boresight_pixels[0,i]+1, boresight_pixels[1,i]-1:boresight_pixels[1,i]+1 ] = max(test_frame_cd) ;hack really need to hard code this and fix the color bar
  
        Range       = [-fov / (3600.*2.),  fov / (3600.*2.)]
        axis_format = { xticks:6, yticks:6, XRange:Range, xtitle:'Degrees', YRange:Range, charsize:0.75}
  
        number_of_colors = 8                                                     ;number of colors / labels on the color bar, loadgood has 8 
  
        P = pos[*,0]+[0, -.05, 0, .1] 
        cgImage, bytscl(alog10(frame), 0., max(alog10(test_frame))), Position = P, AXKEYWORDS=axis_format, /Axes, /KEEP_ASPECT
        cgColorbar, /YLOG, Range=[10.^0., 10.^max(alog10(test_frame))], /Vertical, Position = [.45, P[1], .47, P[3]], $
           Title = Strcompress('Simulated Na Brightness [KR]'), TLocation='Left', /right, divisions = 8, charsize = 0.75, $
           ticknames = strcompress(string(10. ^ (max(alog10(test_frame), /Nan) * findgen(number_of_colors + 1.) / number_of_colors), format = '(F10.1)'), /remove_all)
           
        P = pos[*,1]+[0.015, -.05, 0., .1]
        cgImage, bytscl(alog10(frame_CD), 7., max(alog10(test_frame_CD))), AXKEYWORDS=axis_format, /Axes, /KEEP_ASPECT, Position = P, /noerase, ctindex = 0


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
            
            cgcontour, ob.lat, /onimage, levels = float(STRMID(surface_distribution_this_run, STRPOS(surface_distribution_this_run, ',')+1, STRPOS(surface_distribution_this_run, ']')-STRPOS(surface_distribution_this_run, ',')-1)), $
              LABEL = 0, color = 'red', THICK = .5
            ;cgcontour, ob.lon, /onimage, levels = float(STRMID(Surface_distribution_this_run, STRPOS(Surface_distribution_this_run, '[')+1, STRPOS(Surface_distribution_this_run, ',')-STRPOS(Surface_distribution_this_run, '[')-1)), $
            ;  LABEL = 0, color = 'red', THICK = .5
            ;cgcontour, ob.lon, /onimage, levels = sub_solar_lon*!radeg, color = 'red'
            
            
        cgColorbar, /YLOG, Range=[10.^7., 10.^max(alog10(test_frame_CD), /Nan)], /Vertical, Position = [.935, P[1], .955, P[3]], $
          Title = Strcompress('Column Density [log!D10!N atoms / cm!U-2!N]'), TLocation='Left', /right, divisions = number_of_colors, charsize = 0.75, $
          ticknames = strcompress(string((7. + (max(alog10(test_frame_CD))-7.) * findgen(number_of_colors + 1.) / number_of_colors), format = '(F10.1)'), /remove_all), $
          /noerase, ctindex = 0
        
        
        ;cgColorbar, /Vertical, Position = [.935, P[1], .955, P[3]], Range=[MinValue_CD, MaxValue_CD], $
        ;  Title = Strcompress('Column Density [x 10!U4!N atoms / cm!U-2!N]'), TLocation='Left', /right, $
        ;  ticknames = colorbar_label, divisions = number_of_colors, ticklen = 1, charsize = 0.75, /noerase, ctindex = 0

        P = [pos[0,2]+0.065, pos[1,2], pos[2,3], pos[3,2]]
        cgplot, PDS_JD_Na, DDR_Na, color = 'grey', psym=16, /xstyle, ytitle = 'Na D Brightness [KR]', XTICKUNITS = ['Time'], $
          XTICKFORMAT='LABEL_DATE', xtitle = strmid(utc_string[0], 0, 10)+' UTC', Position=P, yr = [-1., 18.], charsize = 1, /noerase
        cgplot, PDS_JD_Mg, DDR_Mg*10., color = 'red', psym=16, /overplot
          
        cgplot, JD, Tim_RADIANCE_KR, color = 'Black', psym=16, /overplot
        cgplot, JD[0:i], Simulated_UVVS_brightness[0:i], color = 'Blue', psym=16, /overplot   
        cgplot, [JD[i], JD[i]], [0, 20], /overplot
        cgtext, JD[87],15, 'UVVS (PDS DDR)', color = 'grey', charsize = 1
        cgtext, JD[87],13, 'UVVS (Tim''s File)', charsize = 1
        cgtext, JD[87],11, 'Simulation: ' + strmid(Meteor_impact_UTC, 10, 6) + $
                           strcompress( round(brightness_multiplier * sxpar(header, 'UPWARD_F') * 22.989769*1.e-3 / 6.02214e23)) + ' kg', color = 'blue', charsize = 1
        cgtext, JD[87] + .7/24., 9, ' [lon, lat] = ' + strmid(Surface_distribution_this_run, 6), color = 'blue', charsize = 1
        cgtext, JD[87] + .72/24., 7, Plume_Temperature + ' MBF', color = 'blue', charsize = 1
  
      cgPS_Close, width = xs, /PNG;, /Delete_PS ; Convert to PNG file.

      image = Read_PNG(Directory + write_directory + '\movie.png')
      image = image[*,0:xs-1,0:ys-1] ;not sure why this is needed but it is
      ;File_Delete, Directory + write_directory + '\movie.png'
      void = video -> Put(stream, image)
  endfor
  video -> Cleanup
  
stop
end