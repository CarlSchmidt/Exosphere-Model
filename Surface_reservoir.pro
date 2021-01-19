Pro Surface_reservoir

      Plume_Temperature             = '1200K'                    ; temperature of the impact vapour
      Surface_distribution          = 'Global'                    
      loop_times                    = 6.                         ; Bear minimum for any reasonable S/N
      Na_Lofted                     = 1.e25                      ; seems like a lot
      Mg_Lofted                     = 4.*Na_lofted               ; seems like a lot
      Brightness_multiplier_Na      = 1.6
      Brightness_multiplier_Mg      = 0.6

COMMON Output_shared, Plot_range, Output_Size_In_Pixels, Output_Title, Center_in_frame, viewpoint, FOV, N_ticks, Tickstep, Observatory, Above_Ecliptic, Boresight_Pixel, Aperture_Corners
COMMON Model_shared, Body, Ephemeris_time, Seed, Directory, Particle_data, Line_data, Debug

; define the time window to be plotted  
  plot_times       = ['2011-Jun-12 08:00', '2011-Dec-05 06:00'] ; times to plot (UTC) two mercury years
  cspice_utc2et, plot_times[0], start_ET1
  cspice_utc2et, plot_times[1], stop_ET1
  t_intervals = 100 
  ETs = (stop_et1 - start_et1) * findgen(t_intervals)/t_intervals
  
  plot_times       = ['2011-Sep-08 07:00', '2012-Mar-02 06:00'] ; times to plot (UTC) two mercury years
  cspice_utc2et, plot_times[0], start_ET2
  cspice_utc2et, plot_times[1], stop_ET2
  cspice_et2UTC, ets+start_et1, 'C', 0, utcstr1
  cspice_et2UTC, ets+start_et2, 'C', 0, utcstr2
  utcstr = [utcstr1, utcstr2]

  cd, directory
  write_directory             = 'Surface_reservoir\'
  FILE_MKDIR, write_directory, /NOEXPAND_PATH

;; Run the Na simulations at each  time
;  for i = 0, n_elements(utcstr)-1 do begin
;    Generic_Model, Output_title_this_run = write_directory + '\Na_reservoir_'+strcompress(string(i), /remove_all), $
;                   test_particle_this_run = 'Na', Line_this_run = 'Na-D', UTC_this_run = utcstr[i], $
;                   Surface_distribution_this_run = Surface_distribution, Speed_distribution_this_run = 'MBF_'+plume_temperature, $
;                   loop_times_this_run = loop_times, Upward_flux_at_exobase_this_run = Na_Lofted          
;  endfor

  ;-----------------------MOVIE-----------------------------------------
  mpgFilename = Directory+Write_directory+'\Surface_Reservoir_Sequence.mp4'
  framerate   = 10 ; FPS
  
  Reimpact_files = FILE_search('C:\IDL\Generic Model V2\read_write\'+write_directory+'\Na_reservoir_*_Reimpacting_flux_Map.fit', Count = n_files)
  Release_files = FILE_search('C:\IDL\Generic Model V2\read_write\'+write_directory+'\Na_reservoir_*_Release_flux_Map.fit', Count = n_files)
  Reimpacting_flux_Maps = fltarr(360, 180, n_files)
  Release_flux_Maps = fltarr(360, 180, n_files)

  for i = 0, n_files-1 do begin
    Release_flux_Maps[*,*,i]     = mrdfits('C:\IDL\Generic Model V2\read_write\'+write_directory+'\Na_reservoir_'+strcompress(string(i), /remove_all)+'_Release_flux_Map.fit', /silent)
    Reimpacting_flux_Maps[*,*,i] = mrdfits('C:\IDL\Generic Model V2\read_write\'+write_directory+'\Na_reservoir_'+strcompress(string(i), /remove_all)+'_Reimpacting_flux_Map.fit', /silent)
  endfor
  Release_flux_Map     = mean(Release_flux_Maps, dim = 3)
  Reimpacting_flux_Map = mean(Reimpacting_flux_Maps, dim = 3) 

  ;UTC_view_time = '2011-Apr-18 13:00'
  ;UTC_view_time = '2011-Apr-06 18:10' ;TAA 110
  UTC_view_time = '2011-Mar-21 03:55'                ; Mercury 30 deg TAA use for Merkel proposal
  ;UTC_view_time = '2011-Sep-12 18:40:00'                ; Mercury 28 deg TAA use for Merkel proposal
  cspice_str2et, UTC_view_time, ephemeris_time          ; convert UTC to ephemeris time (expressed as the number of ephemeris seconds past J2000)
  
  Display_Surface_Reimpacts, Release_flux_map, Reimpacting_flux_Map, body, ephemeris_time, viewpoint

  ; write the "Reimpacting_Flux_Map" and use it to weight the particles fractional content
    mwrfits, Reimpacting_flux_Map, 'C:\IDL\Generic Model V2\read_write\'+write_directory+'Mean_Reimpacting_Flux_Map.fit', /create 

    Generic_Model, Output_title_this_run = write_directory + 'Cold_Pole_Emission_Day', $
                   test_particle_this_run = 'Na', Line_this_run = 'Na-D', UTC_this_run = UTC_view_time, $
                   Surface_distribution_this_run = 'From_Map', Speed_distribution_this_run = 'MBF_1200K', $
                   loop_times_this_run = 600, Upward_flux_at_exobase_this_run = Na_Lofted, Time_range_this_run = [0.0, 0.25]
stop

  ;--------make the movie--------------------
  mpgFilename = Directory+Write_directory+'\Reimpact_Sequence.mp4'
  framerate   = 10 ; FPS
  loadct, 3
  video      = IDLffVideoWrite(mpgFilename, format='mp4')
  xs = 800
  ys = 600
  stream     = video.AddVideoStream(xs, ys, framerate)

  axis_keywords = {Xrange:[0, 360], Xtitle:'East Longitude', Yrange:[-90, 90], Ytitle:'Latitude', xtickinterval:30, ytickinterval:30}
  for i = 0, n_files-1 do begin
    cgPS_Open, Directory + Write_directory + '\movie.eps', /ENCAPSULATED, /Quiet
    cgDisplay, xs = xs, ys = ys
    device, SET_FONT = 'Helvetica Bold', /TT_FONT

    Release_flux_Maps[*,*,i]     = mrdfits('C:\IDL\Generic Model V2\read_write\'+write_directory+'\Na_reservoir_'+strcompress(string(i), /remove_all)+'_Release_flux_Map.fit')
    Reimpacting_flux_Maps[*,*,i] = mrdfits('C:\IDL\Generic Model V2\read_write\'+write_directory+'\Na_reservoir_'+strcompress(string(i), /remove_all)+'_Reimpacting_flux_Map.fit')

    frame = rebin(Reimpacting_flux_Maps[*,*,i], 180, 90)
    cgimage, frame, minval = 5.e6, maxval = 2.5e7, /axes, /keep, AXKEYWORDS = axis_keywords;, title = 'Atoms deposited cm!U-2!N s!U-1!N'
    cgcolorbar, range = [5.e6, 2.5e7], title = 'Na atoms deposited cm!U-2!N s!U-1!N'
    cgPS_close, width = xs, /PNG;, /Delete_PS ; Convert to PNG file.

    image = Read_PNG(Directory + write_directory + '\movie.png')
    image = image[*,0:xs-1,0:ys-1] ;not sure why this is needed but it is
    void = video -> Put(stream, image)
  endfor
  video -> Cleanup
  
  stop
  ;-----------------movie is made-----------------------------
;    
;    test_frame_Na = mrdfits(Na_files[0], 0, header, /silent) * brightness_multiplier_Na / 1.e3 
;    test_frame_mg = mrdfits(Mg_files[0], 0, header, /silent) * brightness_multiplier_Mg 
;    dimensions    = size(test_frame_Na, /dimensions)
;    xs            = dimensions[0]*8
;    ys            = dimensions[1]*6
; 
;  ; Format layout and the time axis
;    pos        = cgLayout([2,2], xgap = 0, ixmargin = 0, iymargin = 1, oymargin = 5, oxmargin = 0) 
;    date_label = LABEL_DATE(DATE_FORMAT = ['%H:%I'])
;    video      = IDLffVideoWrite(mpgFilename, format='mp4')
;    stream     = video.AddVideoStream(xs, ys, framerate)
;    !P.font    = 1  
;    Simulated_UVVS_brightness_Na = MAKE_ARRAY(N_elements(DDR_Na), /float, VALUE = !values.f_nan) 
;    Simulated_UVVS_brightness_Mg = MAKE_ARRAY(N_elements(DDR_Mg), /float, VALUE = !values.f_nan) 
;    loadgood
;
;  for i = 0, N_elements(DDR_Na)-1 do begin
;    
;    ; Find the nearest Mg Frame
;      junk = min(abs(PDS_ET_Na[i] - PDS_ET_Mg[Key_frames_Mg]), which_Mg)
;      nearest_Mg  = Key_frames_Mg[which_Mg]
;  
;    ; only plot "key frames"
;      junk = where(key_frames_Na eq i, count, /NULL)
;      if count eq 0 then continue
;    
;    ; And only plate files that exist since some measurements are pre-impact  
;      if not (FILE_TEST('C:\IDL\Generic Model V2\read_write\'+write_directory+'\Na_Meteor_frame_'+strcompress(string(i), /remove_all)+'.fit') or $
;              FILE_TEST('C:\IDL\Generic Model V2\read_write\'+write_directory+'\Mg_Meteor_frame_'+strcompress(string(i), /remove_all)+'.fit')) then continue
;
;    ; We are going to create high-resolution PNG images, which we will add to the video stream.
;    ; We create the high-resolution PNG image from PostScript intermediate files.
;      cgPS_Open, Directory + Write_directory + '\movie.eps', /ENCAPSULATED, /Quiet
;      cgDisplay, xs = xs, ys = ys
;      device, SET_FONT = 'Helvetica Bold', /TT_FONT
;
;      ; read the model images
;        Meteor_frame_Na = mrdfits(Directory + write_directory + '\Na_Meteor_frame_'+strcompress(string(i), /remove_all)+'.fit', 0, Na_header, /silent) 
;        Meteor_frame_Mg = mrdfits(Directory + write_directory + '\Mg_Meteor_frame_'+strcompress(string(Nearest_Mg), /remove_all)+'.fit', 0, Mg_Header, /silent)
;        platescale      = float(sxpar(Na_header, 'FOV')) / float(sxpar(Na_header, 'NAXIS1'))  ; "/pix image platescale
;        
;        if (N_elements(Meteor_frame_Na) eq 1) or (N_elements(Meteor_frame_Mg) eq 1) then continue
;      
;      ; read the UVVS pointing info  
;        Pointing_Na      = mrdfits(Directory + write_directory + '\Na_Meteor_frame_'+strcompress(string(i), /remove_all)+'.fit', 2, /silent)
;        Pointing_Mg      = mrdfits(Directory + write_directory + '\Mg_Meteor_frame_'+strcompress(string(Nearest_Mg), /remove_all)+'.fit', 2, /silent)
;
;        ;frame       = (background_frame + Meteor_frame) / 1.e3 ; Convert brightness from R to KR
;        frame_Na      =  brightness_multiplier_Na * Meteor_frame_Na / 1.e3 ; Convert brightness from R to KR
;        frame_Mg      =  brightness_multiplier_Mg * Meteor_frame_Mg 
;      
;      ; Log the simulated brightness, assumed to be the mean of the brightness at the apertures 4 corner positions 
;        ;Simulated_UVVS_brightness_Na[i] = mean( frame_Na[ boresight_pixels[0]-1:boresight_pixels[0]+1, boresight_pixels[1]-1:boresight_pixels[1]+1 ] )
;        Simulated_UVVS_brightness_Na[i] = mean(interpolate(frame_Na, Pointing_Na.APERTURE_CORNERS[0,*], Pointing_Na.APERTURE_CORNERS[1,*])) 
;        Simulated_UVVS_brightness_Mg[Nearest_Mg] = mean(interpolate(frame_Mg, Pointing_Mg.APERTURE_CORNERS[0,*], Pointing_Mg.APERTURE_CORNERS[1,*])) 
;  
;        Range       = [-fov / (3600.*2.),  fov / (3600.*2.)]
;        axis_format = { xticks:6, yticks:6, XRange:Range, xtitle:'Degrees', YRange:Range, charsize:0.75}
;        number_of_colors = 8                                                     ;number of colors / labels on the color bar, loadgood has 8 
;  
;        P = pos[*,0]+[0, -.03, 0, .08] 
;        cgImage, bytscl(alog10(frame_Na), 0., max(alog10(test_frame_Na))), Position = P, AXKEYWORDS=axis_format, /Axes, /KEEP_ASPECT
;        cgColorbar, Range=[0., max(alog10(test_frame_Na))], /Vertical, Position = [.445, P[1], .465, P[3]], $
;          Title = Strcompress('Simulated Na Brightness Log!D10!N[kR]'), TLocation='Left', /right, divisions = number_of_colors, charsize = 0.75
;        
;        ; Mark the UVVS aperture projection, the coordinates here need to be in degrees  
;          cgPolygon, [transpose(Pointing_Na.APERTURE_CORNERS[0,*]), Pointing_Na.APERTURE_CORNERS[0,0]]*platescale/3600.+range[0], $
;                     [transpose(Pointing_Na.APERTURE_CORNERS[1,*]), Pointing_Na.APERTURE_CORNERS[1,0]]*platescale/3600.+range[0], COLOR='Snow', /fill
;
;         ;----------------------------Generate the Latitude & Longitude grid to overlay on the column density----------------------------------
;         ; Find the angular size
;           CSPICE_SPKEZR, body, PDS_ET_Na[i], 'J2000', 'LT', viewpoint, BODY_state, ltime
;           R_M = 206264.806 * atan(radii[0] / norm(BODY_state[0:2]))          ; Radius of the body in arcsec
;
;         ; Find the sub-observer planetographic longitude and latitude
;           cspice_subpnt, 'Near point: ellipsoid', Body, PDS_ET_Na[i], 'IAU_'+body, 'LT+S', viewpoint, $
;             spoint, trgepc, srfvec
;           f = ( radii[0]-radii[2] ) / radii[0]                               ; flatness parameter, mercury is round so this is zero
;           cspice_recpgr, body, spoint, radii[0], f, spglon, spglat, spgalt
;           lat_se = spglat * cspice_dpr()                                     ; planetographic latitude of the sub-observer point on the surface. (0 at equator, 90 at north pole, -90 at south pole)
;           lon_se = spglon * cspice_dpr()                                     ; planetographic longitude of the sub-observer point on the surface. (0 to 360 west longitude).
;           ;if keyword_set(debug) then print, 'Sub-observer Lat & Lon  =', lat_se, lon_se, ' Angular diameter of ', body, ' =', 2.*R_M, ' arcsec'
;
;         ; setup the x-y grid image
;           xdim     = float(sxpar(Na_header, 'NAXIS1'))
;           ydim     = float(sxpar(Na_header, 'NAXIS2'))
;           ctr_xpix = xdim / 2.
;           ctr_ypix = ydim / 2.
;           pix2km   = platescale * (1./R_M) * radii[0]                        ; "/pix * radii/" * km/radii
;           x        = (dindgen(xdim) - ctr_xpix) * pix2km
;           y        = (dindgen(xdim) - ctr_ypix) * pix2km
;           xsq      = dblarr(xdim, ydim)                                      ; xsq is an image the size of the calib-img where the values are the horizontal distance from the "center" in km
;           ysq      = dblarr(xdim, ydim)                                      ; ysq is an image the size of the calib-img where the values are the vertical distance from the center.
;           for j = 0, xdim-1 do begin
;             xsq[*,j] = x
;             ysq[j,*] = y
;           endfor
;
;         ; Get the planetographic lon/lat at each pixel.
;         ; This grid looks best if make the arrays big, rotate things to our "Sunward = -X, Observer-Planet = +Z" geometry with interpolation, then shrink shrink them again.
;           ysq = rebin(ysq, 1024, 1024)
;           xsq = rebin(xsq, 1024, 1024)
;           ob     = deprob(xsq, ysq, radii[0], radii[2], lat_se, lon_se)
;           ob.lon[where(ob.lon eq -666, /Null)] = !values.F_nan
;           ob.lat[where(ob.lat eq -666, /Null)] = !values.F_nan
;           ob.lon = rot(ob.lon, sub_solar_lon*!radeg - lon_se - 90., /interp)
;           ob.lat = rot(ob.lat, sub_solar_lon*!radeg - lon_se - 90., /interp)
;           ob.lon = rebin(ob.lon, 1024, 1024)
;           ob.lat = rebin(ob.lat, 1024, 1024)
;
;         ; Mark the impacte site's lat and lon
;           cgcontour, ob.lat, /onimage, levels = float(STRMID(Surface_distribution, STRPOS(Surface_distribution, ',')+1, STRPOS(Surface_distribution, ']')-STRPOS(Surface_distribution, ',')-1)), $
;             color = 'red', THICK = 1.5, LABEL = 0
;           cgcontour, ob.lon, /onimage, levels = float(STRMID(Surface_distribution, STRPOS(Surface_distribution, '[')+1, STRPOS(Surface_distribution, ',')-STRPOS(Surface_distribution, '[')-1)), $
;             color = 'red', THICK = 1.5, LABEL = 0, OLEVELS = OLEVELS
;           ;cgcontour, ob.lon, /onimage, levels = sub_solar_lon*!radeg, color = 'YELLOW', THICK = 0.5 ; Mark the Sub-Solar Point
;
;         ; Mark the lat & lon grid
;           lat_contours = indgen(17)*10 - 80                                  ; every 10 deg lat
;           lon_contours = indgen(12)*30                                       ; every 30 deg lon
;           cgcontour, ob.lat, /onimage, levels = lat_contours, LABEL = 0, color = 'Snow', THICK = .5
;           cgcontour, ob.lon, /onimage, levels = lon_contours, LABEL = 1, color = 'Snow', THICK = .5
;           
;        P = pos[*,1]+[0.015, -.03, 0., .08]
;        cgImage, bytscl(alog10(frame_Mg), 2., max(alog10(test_frame_Mg))), AXKEYWORDS=axis_format, /Axes, /KEEP_ASPECT, Position = P, /noerase
;        cgColorbar, Range=[2., max(alog10(test_frame_Mg), /Nan)], /Vertical, Position = [.93, P[1], .95, P[3]], /noerase, $
;                    Title = Strcompress('Simulated Mg Brightness Log!D10!N[R]'), TLocation='Left', /right, divisions = number_of_colors, charsize = 0.75
;
;        ; Mark the UVVS aperture projection, the coordinates here need to be in degrees
;          cgPolygon, [transpose(Pointing_Mg.APERTURE_CORNERS[0,*]), Pointing_Mg.APERTURE_CORNERS[0,0]]*platescale/3600.+range[0], $
;                     [transpose(Pointing_Mg.APERTURE_CORNERS[1,*]), Pointing_Mg.APERTURE_CORNERS[1,0]]*platescale/3600.+range[0], COLOR='Snow', /fill
;                     
;            ;----------------------------Generate the Latitude & Longitude grid to overlay on the column density----------------------------------
;            ; Find the angular size
;              CSPICE_SPKEZR, body, PDS_ET_Mg[nearest_Mg], 'J2000', 'LT', viewpoint, BODY_state, ltime
;              R_M = 206264.806 * atan(radii[0] / norm(BODY_state[0:2]))          ; Radius of the body in arcsec
;
;            ; Find the sub-observer planetographic longitude and latitude
;              cspice_subpnt, 'Near point: ellipsoid', Body, PDS_ET_Mg[nearest_Mg], 'IAU_'+body, 'LT+S', viewpoint, $
;                spoint, trgepc, srfvec
;              f = ( radii[0]-radii[2] ) / radii[0]                               ; flatness parameter, mercury is round so this is zero
;              cspice_recpgr, body, spoint, radii[0], f, spglon, spglat, spgalt
;              lat_se = spglat * cspice_dpr()                                     ; planetographic latitude of the sub-observer point on the surface. (0 at equator, 90 at north pole, -90 at south pole)
;              lon_se = spglon * cspice_dpr()                                     ; planetographic longitude of the sub-observer point on the surface. (0 to 360 west longitude).
;              ;if keyword_set(debug) then print, 'Sub-observer Lat & Lon  =', lat_se, lon_se, ' Angular diameter of ', body, ' =', 2.*R_M, ' arcsec'
;    
;            ; setup the x-y grid image
;              xdim     = float(sxpar(Mg_header, 'NAXIS1'))
;              ydim     = float(sxpar(Mg_header, 'NAXIS2'))
;              ctr_xpix = xdim / 2.
;              ctr_ypix = ydim / 2.
;              pix2km   = platescale * (1./R_M) * radii[0]                        ; "/pix * radii/" * km/radii
;              x        = (dindgen(xdim) - ctr_xpix) * pix2km
;              y        = (dindgen(ydim) - ctr_ypix) * pix2km
;              xsq      = dblarr(xdim, ydim)                                      ; xsq is an image the size of the calib-img where the values are the horizontal distance from the "center" in km
;              ysq      = dblarr(xdim, ydim)                                      ; ysq is an image the size of the calib-img where the values are the vertical distance from the center.
;              for j = 0, xdim-1 do begin
;                xsq[*,j] = x
;                ysq[j,*] = y
;              endfor
;    
;            ; Get the planetographic lon/lat at each pixel. 
;            ; This grid looks best if make the arrays big, rotate things to our "Sunward = -X, Observer-Planet = +Z" geometry with interpolation, then shrink shrink them again.
;              ysq = rebin(ysq, 1024, 1024)
;              xsq = rebin(xsq, 1024, 1024)
;              ob     = deprob(xsq, ysq, radii[0], radii[2], lat_se, lon_se)
;              ob.lon[where(ob.lon eq -666, /Null)] = !values.F_nan
;              ob.lat[where(ob.lat eq -666, /Null)] = !values.F_nan             
;              ob.lon = rot(ob.lon, sub_solar_lon*!radeg - lon_se - 90., /interp)
;              ob.lat = rot(ob.lat, sub_solar_lon*!radeg - lon_se - 90., /interp)
;              ob.lon = rebin(ob.lon, 1024, 1024)
;              ob.lat = rebin(ob.lat, 1024, 1024)
;      
;            ; Mark the impacte site's lat and lon
;              cgcontour, ob.lat, /onimage, levels = float(STRMID(Surface_distribution, STRPOS(Surface_distribution, ',')+1, STRPOS(Surface_distribution, ']')-STRPOS(Surface_distribution, ',')-1)), $
;                color = 'red', THICK = 1.5, LABEL = 0
;              cgcontour, ob.lon, /onimage, levels = float(STRMID(Surface_distribution, STRPOS(Surface_distribution, '[')+1, STRPOS(Surface_distribution, ',')-STRPOS(Surface_distribution, '[')-1)), $
;                color = 'red', THICK = 1.5, LABEL = 0, OLEVELS = OLEVELS
;              ;cgcontour, ob.lon, /onimage, levels = sub_solar_lon*!radeg, color = 'YELLOW', THICK = 0.5 ; Mark the Sub-Solar Point
;            
;            ; Mark the lat & lon grid  
;              lat_contours = indgen(17)*10 - 80                                  ; every 10 deg lat
;              lon_contours = indgen(12)*30                                       ; every 30 deg lon
;              cgcontour, ob.lat, /onimage, levels = lat_contours, LABEL = 0, color = 'Snow', THICK = .5
;              cgcontour, ob.lon, /onimage, levels = lon_contours, LABEL = 1, color = 'Snow', THICK = .5
;
;        P = [pos[0,2]+0.065, pos[1,2], pos[2,3]-0.065, pos[3,2]]
;        cgplot, PDS_JD_Na, DDR_Na, color = 'Black', psym=16, /xstyle, YSTYLE=9, ytitle = 'Na D Brightness [kR]', XTICKUNITS = ['Time'], $
;          XTICKFORMAT='LABEL_DATE', xtitle = strmid(Na_UTC_string[0], 0, 10)+' UTC', Position=P, yr = [-1., 18.], charsize = 1, /noerase
;        cgplot, PDS_JD_Mg, DDR_Mg*10., color = 'red', psym=16, XRANGE = !X.CRANGE, Position=P, charsize = 1, /noerase, ystyle = 16, XTICKUNITS = ['Time'], $
;          XTICKFORMAT='LABEL_DATE', yr = [-1., 18.]
;        cgAXIS, YAXIS=1, YRANGE = !Y.CRANGE*100., YSTYLE = 1, YTITLE = 'Mg 2853A Brightness [R]', charsize = 1., color = 'red'
;
;        cgplot, PDS_JD_Na[0:i], Simulated_UVVS_brightness_Na[0:i], color = 'grey', psym=14, /overplot;, symsize = 0.75      
;        cgplot, PDS_JD_Mg[0:Nearest_Mg], Simulated_UVVS_brightness_Mg[0:Nearest_Mg]/100., color = 'pink', psym=14, /overplot;, symsize = 0.75   
;
;        cgplot, [PDS_JD_Na[i], PDS_JD_Na[i]], [0, 20], /overplot
;        cgtext, PDS_JD_Na[87],15, 'UVVS Na (SNR > 3)', color = 'Black', charsize = 1
;        cgtext, PDS_JD_Na[87],13, 'UVVS Mg (SNR > 3)', Color = 'red', charsize = 1
;        cgtext, PDS_JD_Na[87],11, 'Simulation: ' + strmid(Meteor_impact_UTC, 10, 6) + ' [lon, lat] = ' + strmid(Surface_distribution, 6), color = 'black', charsize = 1
;        cgtext, PDS_JD_Na[87] + .72/24., 9, Plume_Temperature + ' MBF', color = 'black', charsize = 1
;        cgtext, PDS_JD_Na[87] + .70/24., 7, strcompress( round(brightness_multiplier_Na * sxpar(Na_header, 'UPWARD_F') * 22.989769*1.e-3 / 6.02214e23)) + ' kg (Na)', color = 'grey', charsize = 1
;        cgtext, PDS_JD_Na[87] + 1.5/24., 7, strcompress( round(brightness_multiplier_Mg * sxpar(Mg_header, 'UPWARD_F') * 24.305*1.e-3 / 6.02214e23)) + ' kg (Mg)', color = 'pink', charsize = 1
;
;      cgPS_Close, width = xs, /PNG;, /Delete_PS ; Convert to PNG file.
;      image = Read_PNG(Directory + write_directory + '\movie.png')
;      image = image[*,0:xs-1,0:ys-1] ;not sure why this is needed but it is
;      void = video -> Put(stream, image)
;  endfor
;  video -> Cleanup
;
;  ; Compare the results by correlation
;    correl_Na = correlate(Simulated_UVVS_brightness_Na[key_frames_Na], DDR_Na[key_frames_Na])
;    correl_Mg = correlate(Simulated_UVVS_brightness_Mg[key_frames_Mg], DDR_Mg[key_frames_Mg])
;    Print, 'correlation of results Na & Mg [Key_frames]:', correl_Na, correl_Mg
;    save, correl_Na, correl_Mg, filename = Directory+ 'Correlation_results\' + write_directory + '_correlation_Hi-res.sav'
end