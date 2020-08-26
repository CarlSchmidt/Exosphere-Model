pro output_display, loc, g, atoms_per_packet, Image_type, loop_number, part, Label_North = Label_North, Label_phase = Label_phase, $
    Label_time = Label_time 

; REVISION HISTORY:
;    Written:                     C. Schmidt, July 2009
;    -Vectorized some, added 'Above Ecliptic' & 'Moon Spot' viewpoints
;    -Changed Plot_Range to FOV: from planetary radii to arcseconds ---> need to figure out above ecliptic scale in this case!
;                                 C. Schmidt, Nov 2018

COMMON Output_shared, Plot_range, Output_Size_In_Pixels, Output_Title, Center_in_frame, viewpoint
COMMON Model_shared, Body, Ephemeris_time, Seed, Directory, Particle_data, Line_data, Debug


FOV = 60. ;Field of view in arcseconds FIX FOR NOW HACK HACK!!!
N_ticks = 10.
ang_platescale = float(Output_Size_In_Pixels[0]) / float(FOV) ;plate scale of the output image in PIXELS per arcsec 
Observatory = 'McDonald' 

if part eq 1 then begin

  ; Looking down over the ecliptic is a special case
    if viewpoint eq 'Above Ecliptic' then begin
        cspice_pxform, "J2000", "GSE", ephemeris_time, J2000_to_GSE_transform_matrix ;'Z' axis is normal to the ecliptic now.
        loc_prime = TRANSPOSE(J2000_to_GSE_transform_matrix) # loc[0:2,*]
    endif
  
  ; Looking down over the ecliptic is a special case
    if viewpoint eq 'Moon Spot' then begin
        ; Get the body-fixed non-inertial IAU_Earth Coordinates of the Observatory in Cartesian
          cspice_bodvrd, 'EARTH', 'RADII', 3, radii
          flat = (radii[0] - radii[2])/radii[0]
          OBSERVATORY_cs, 'McDonald', obs_struct ;longitude is in degrees *West*, need to switch to East Longit for J2000 conversion via cspice_georec
          cspice_georec, (360. - obs_struct.longitude) * CSPICE_RPD(), obs_struct.latitude * CSPICE_RPD(), obs_struct.altitude / 1.e3, radii[0], flat, obs_IAU_Earth
      
        ; convert to J2000
          cspice_pxform, 'IAU_EARTH', 'J2000', ephemeris_time, Earth_to_J2000_transform_matrix 
          obs_J2000 = TRANSPOSE(Earth_to_J2000_transform_matrix) # obs_IAU_Earth

        ; Get the RA and Dec of the anti-sunward vector Sun 
          cspice_spkpos, 'Sun', ephemeris_time, 'J2000', 'LT', 'Earth', Sun_Earth_vector, light_time 
          cspice_recrad, -(Sun_Earth_vector-obs_J2000), radius, anti_Sun_RA, anti_Sun_dec
          cspice_recrad, (Sun_Earth_vector-obs_J2000), radius, Sun_RA, Sun_dec
          ;print, RA*!radeg, Dec*!radeg ; ---> sunward vector verified within 0.055" of the J2000 Horizons result from McDonald Obs
          
        ; Now get the RA and Dec of all the particles
          ; Test particles are in the J2000 frame with body/Lunar center as the origin 
          ; Move Loc array origin to observatory, still in J2000
          cspice_spkpos, body, ephemeris_time, 'J2000', 'LT', 'Earth', Body_Earth_vector, light_time        
          Loc_WRT_obs_J2000 = loc[0:2,*] - rebin(Body_Earth_vector, 3, N_elements(loc[0,*])) - rebin(obs_J2000, 3, N_elements(loc[0,*])) 
          cspice_recrad, Loc_WRT_obs_J2000, Loc_radius_WRT_obs, Loc_RA_WRT_obs, Loc_dec_WRT_obs
          
          ;looking sunward
          window, 0
          x = (Loc_RA_WRT_obs - Sun_RA)*CSPICE_RPD()
          Y = (Loc_dec_WRT_obs - Sun_dec)*CSPICE_RPD()
          cgplot, x, y, xr = [-180,180], yr = [-90,90], psym = 4    
          
          ;looking anti-sunward
          window, 1 
          x = (Loc_RA_WRT_obs - anti_Sun_RA)*CSPICE_RPD()
          Y = (Loc_dec_WRT_obs - anti_Sun_dec)*CSPICE_RPD()
          cgplot, x, y, xr = [-40,40], yr = [-40,40], psym = 4         
    stop
    endif      
    
  ; Normal Cases: One object viewed from the center of another object
    if (viewpoint ne 'Above Ecliptic' and viewpoint ne 'Moon Spot') then begin
      ; Generate a coordinate system where the Z axis lies along the Observer-Body line, and the *projected* anti-sunward tail is along the -x axis.
        cspice_spkpos, 'Sun', ephemeris_time, 'J2000', 'NONE', body, Sun_Body_vector, light_time 
        cspice_spkpos, viewpoint, ephemeris_time, 'J2000', 'NONE', body, Observer_Body_vector, light_time 
    
      ; Now that we have the Observer to Body and Sun to Body vectors, we can define a coordinate system in which:
      ; The +Z axis is the Observer to Body vector and The Sun-Body vector lies in the XZ plane.
      
      ; Finding the rotation matrix to transform particles from the current J2000, Solar System Barycenter centered plane
      ; into the desired coordinate frame.
        cspice_twovec, Observer_Body_vector, 3, Sun_Body_vector, 1, alignment_matrix  ;defines -Observer_Body_vector as 'Z' axis
      
      ; Put the packet positions into this new "aligned" coordinate frame: loc_prime
        loc_prime = TRANSPOSE(alignment_matrix) # loc[0:2,*]
    endif  
    
  ; Find the body's radius in Km
    cspice_bodvrd, Body, 'RADII', 3, Body_radius 
    Body_radius = Body_radius[0]

;  ; Recall the Northward unit vector and transform it from the J2000 frame to the "tail aligned frame"
;    restore, strcompress(directory+'North_J2000_frame.sav'), DESCRIPTION = North_vector_info   
;    cspice_mxv, alignment_matrix, North_J2000_frame, North_tail_aligned_frame
;    save, North_tail_aligned_frame, filename = strcompress(directory+'North_tail_aligned_frame.sav'), DESCRIPTION = North_vector_info  
    
  ; Place the body center in pixel space
    x_center = center_in_frame[0]*Output_Size_In_Pixels[0] 
    y_center = center_in_frame[1]*Output_Size_In_Pixels[1] 
    z_center = center_in_frame[2]*Output_Size_In_Pixels[2]      
    
  ;------------------------------------------------BUILD AN IMAGE VEIWED FROM A FINITE RADIUS USING AN ANGULAR PIXEL SIZE-----------------------------------------------------------------
     Model_Image_R  = fltarr(Output_Size_In_Pixels[0], Output_Size_In_Pixels[1]) ;Image in Rayleighs in the XY plane to be written 
     Model_Image_CD = fltarr(Output_Size_In_Pixels[0], Output_Size_In_Pixels[1]) ;Image in cgs column density in the XY plane to be written 
     ; Get the RA and Dec of the 'Body' and all the particles from the observers location...       
        ; Get the body-fixed non-inertial IAU_Earth Coordinates of the Observatory in Cartesian
          cspice_bodvrd, 'EARTH', 'RADII', 3, radii
          flat = (radii[0] - radii[2])/radii[0]
          OBSERVATORY_cs, observatory, obs_struct ;longitude is in degrees *West*, need to switch to East Longit for J2000 conversion via cspice_georec
          cspice_georec, (360. - obs_struct.longitude) * CSPICE_RPD(), obs_struct.latitude * CSPICE_RPD(), obs_struct.altitude / 1.e3, radii[0], flat, obs_IAU_Earth
      
        ; convert to J2000
          cspice_pxform, 'IAU_EARTH', 'J2000', ephemeris_time, Earth_to_J2000_transform_matrix 
          obs_J2000 = TRANSPOSE(Earth_to_J2000_transform_matrix) # obs_IAU_Earth

        ; Get the RA and Dec of the Body 
          cspice_spkpos, Body, ephemeris_time, 'J2000', 'LT', 'Earth', body_Earth_vector, light_time 
          cspice_recrad, (body_Earth_vector - obs_J2000), Body_distance_WRT_obs, Body_RA_WRT_obs, Body_dec_WRT_obs
          ;RA and Dec from an observatory verified accurate within ~0.055" via JPL Horizons (for Sun & Mercury) 
        
        ; Now get the RA and Dec of all the particles
          ; Test particles are in the J2000 frame with body center as the origin 
          cspice_spkpos, body, ephemeris_time, 'J2000', 'LT', 'Earth', Body_Earth_vector, light_time        
          cspice_recrad, loc[0:2,*] + rebin(Body_Earth_vector, 3, N_elements(loc[0,*])) - rebin(obs_J2000, 3, N_elements(loc[0,*])), $
                         Loc_distance_WRT_obs, Loc_RA_WRT_obs, Loc_dec_WRT_obs

    ; From the Delta RA and Delta Dec, build and image...               
          body_ang_radius = atan(Body_radius / Body_distance_WRT_obs)*3600.d/CSPICE_RPD()
          delta_RA        = (Loc_RA_WRT_obs - Body_RA_WRT_obs)*3600.d/CSPICE_RPD()
          delta_Dec       = (Loc_dec_WRT_obs - Body_dec_WRT_obs)*3600.d/CSPICE_RPD()


          points = (2. * !PI / 999.0) * FINDGEN(1000) 
          x_circle = x_center + ang_platescale * COS(points) * body_ang_radius
          y_circle = y_center + ang_platescale * SIN(points) * body_ang_radius
          disc = polyfillv(x_circle, y_circle, Output_Size_In_Pixels[0], Output_Size_In_Pixels[1])
          if disc[0] ne -1 then begin ;if the plate scale is such that a pixel fits on the bodys disc, get the coordinates of the disc.   
            body_resolved = 'true'                                    ; Signifies there will be pixels on the body's disk 
            dims = SIZE(Model_Image_R, /DIMENSIONS) 
            disc_coordinates = array_indices(dims, disc, /DIMENSIONS) ; Coordinates of pixels occuped by body's disc in x and y
          endif else body_resolved = 'false'     
      ;put the packets into the pixel bins
      for i = 0, n_elements(loc[0,*])-1. do begin 
        if (loc[4,i] eq 0.) then continue ;To save computation time, skip particles that have been lost via surface collisions
        x = delta_RA[i] * ang_platescale + x_center
        y = delta_Dec[i] * ang_platescale + y_center 

       ;if the packet location is on the disc show only packets between the body and the observer 
        if body_resolved then here = where(disc_coordinates[0,*] eq fix(x) and disc_coordinates[1,*] eq fix(y), count) else count = 0.
        IF count NE 0 THEN on_disc=1. ELSE on_disc=0. 
        if (on_disc eq 1) and (Loc_distance_WRT_obs[i] gt Body_distance_WRT_obs) then continue ;skip packets behind the disc      
          if x ge 0 and x lt Output_Size_In_Pixels[0] and y ge 0 and y lt Output_Size_In_Pixels[1] then begin ;allocate x,y positions to pixels in the 2d Model_Image_R  if within the plot window 
              ;output an image of the column density
              Model_Image_CD[x,y] = Model_Image_CD[x,y] + ((loc[4,i]) * (atoms_per_packet) * ( ang_platescale*body_ang_radius/(1.e5*body_radius) )^2.)
              ; atoms / cm^2 = (fractional packet content) * (atoms per packet) * ((pixels/cm)^2)  
              ; where pixels/cm = (pixels/arcsec) * (arcsec/cm)
            
              ;output an image in rayleighs
              Model_Image_R[x,y] = Model_Image_R[x,y] + ((loc[4,i]) * (atoms_per_packet) * g[i] * 1.e-6 * ( ang_platescale*body_ang_radius/(1.e5*body_radius) )^2.)              
              ;(fractional packet content) * (atoms per packet) * (photons/(atom*s)) * (10^-6 for Rayleighs) * ((pixels/cm)^2) 
          endif
      endfor
      
;      TICLABELS, Body_dec_WRT_obs*!radeg - FOV/(2.*3600.), N_ticks+1, fov/(60.*(N_ticks)), declabs, delta = 1
;      TICLABELS, Body_RA_WRT_obs*!radeg - (FOV*cos(Body_dec_WRT_obs)/(2.*3600.)), N_ticks+1., (fov*cos(Body_dec_WRT_obs)/(60.*N_ticks)) / 15., RAlabs, /RA, delta = 1
;
;      window, 0, xs = Output_Size_In_Pixels[0]*3, ys = Output_Size_In_Pixels[1]*3
;      axis_format = {XTicks:N_elements(RALABS)-1, XTickname:RALABS, Yticks:N_elements(DecLABS)-1, Ytickname:DecLABS, CHARSIZE:1.2, TICKLEN:-0.01}
;      cgimage, Model_Image_R, AXKEYWORDS=axis_format, /Axes, /KEEP_ASPECT_RATIO

    ;;;--------------------------------------------------end testing--------------------------------------------------------

;    
;  ; Get the coordinates in the x-y plane for the body's disc
;    platescale = float(float(Output_Size_In_Pixels[0]) / (float(plot_range) * body_radius)) ;plate scale of the output image in PIXELS per KM
;    Model_Image_R  = fltarr(Output_Size_In_Pixels[0], Output_Size_In_Pixels[1]) ;Image in Rayleighs in the XY plane to be written 
;    Model_Image_CD = fltarr(Output_Size_In_Pixels[0], Output_Size_In_Pixels[1]) ;Image in cgs column density in the XY plane to be written 
;    Model_Image_V  = fltarr(Output_Size_In_Pixels[0], Output_Size_In_Pixels[1]) ;Image in km/s to be written
;    optically_thick = Model_Image_R  ;use this to generate a seperate image for optically thick pixels, then replace where needed in the image .      
;    points = (2. * !PI / 999.0) * FINDGEN(1000) 
;    x_circle = x_center + platescale * COS(points) * body_radius
;    y_circle = y_center + platescale * SIN(points) * body_radius
;    disc = polyfillv(x_circle, y_circle, Output_Size_In_Pixels[0], Output_Size_In_Pixels[1])
;    if disc[0] ne -1 then begin ;if the plate scale is such that a pixel fits on the bodys disc, get the coordinates of the disc.   
;      body_resolved = 1. 
;      dims = SIZE(Model_Image_R, /DIMENSIONS) 
;      disc_coordinates = array_indices(dims, disc, /DIMENSIONS) ;coordinates of pixels occupes by Body's disc in x and y
;    endif else body_resolved = 0.     
;    ;put the packets into the pixel bins
;    for i = 0, n_elements(loc[0,*])-1. do begin 
;      if (loc[4,i] eq 0.) then continue ;To save computation time, skip particles that have been lost via surface collisions
;      
;      x = loc_prime[0,i] * platescale + x_center ;(distance from Body in km) * (pixels per km) + (Body center pixel location, x) 
;      y = loc_prime[1,i] * platescale + y_center ;(distance from Body in km) * (pixels per km) + (Body center pixel location, y) 
;      z = loc_prime[2,i] * platescale + z_center ;(distance from Body in km) * (pixels per km) + (Body center pixel location, z)      
;;      v_x = loc_prime[5,i] ; km/s
;;      v_y = loc_prime[6,i] ; km/s
;;      v_z = loc_prime[7,i] ; km/s
;  
;      ;if the packet location is on the disc show only packets in the negative z (between the body and the observer) 
;      if body_resolved then here = where(disc_coordinates[0,*] eq fix(x) and disc_coordinates[1,*] eq fix(y), count) else count = 0.
;      IF count NE 0 THEN on_disc=1. ELSE on_disc=0. 
;      if (on_disc eq 1) and (loc_prime[2,i] lt 0.) then continue ;skip packets behind the disc      
;        if x ge 0 and x lt Output_Size_In_Pixels[0] and y ge 0 and y lt Output_Size_In_Pixels[1] then begin ;allocate x,y positions to pixels in the 2d Model_Image_R  if within the plot window 
;            ;output an image of the column density
;            Model_Image_CD[x,y] = Model_Image_CD[x,y] + ((loc[4,i]) * (atoms_per_packet) * (platescale/1.e5)^2.)
;            ;(fractional packet content) * (atoms per packet) * ((pixels/cm)^2)  
;          
;            ;output an image in rayleighs
;            Model_Image_R[x,y] = Model_Image_R[x,y] + ((loc[4,i]) * (atoms_per_packet) * g[i] * 1.e-6 * (platescale/1.e5)^2.)              
;            ;(fractional packet content) * (atoms per packet) * (photons/(atom*s)) * (10^-6 for Rayleighs) * ((pixels/cm)^2) 
;            ;optically_thick = optically_thick[x,y]+((loc(4,i))*(atoms_per_packet)*(platescale/(body_radius*1.e5))^2.) ;column density in atoms/cm^2          
;            ;(fractional packet content) * (atoms per packet) * ((pixels/cm)^2)  
;            ;where optically thick apply a correction:                 
;         endif
;    endfor
;  
;  ;        ;vectorized version, takes even longer!        
;  ;        xs = Output_Size_In_Pixels[0] & ys = Output_Size_In_Pixels[1]
;  ;        Model_Image_R    = FLTARR(xs, ys) ;Model_Image_R  in the XY plane to be written
;  ;        x       = loc_prime[0,*]*platescale + x_center ;(distance from target in Rm) * (pixels per m) + (1/2 the Model_Image_R  length in pixels to center the plot)
;  ;        y       = loc_prime[1,*]*platescale + y_center ;(distance from target in Rm) * (pixels per m) + (1/2 the Model_Image_R  length in pixels to center the plot) 
;  ;        rad     = SQRT((x - x_center)^2 + (y - y_center)^2) ;Radial distance of each particle from body center
;  ;        bh_disc = (rad LT (platescale)) AND (loc_prime[2,*] GT 0.) ;Determine if particles are behind the body disc
;  ;        off_pic = ((x LT 0) OR (x GT xs)) OR ((y LT 0) OR (y GT ys)) ;Determine if particles are beyond the edge of the picture
;  ;        
;  ;        FOR i = 0, n_elements(loc[0,*]) - 1 DO BEGIN
;  ;          IF (bh_disc OR off_pic)[i] THEN CONTINUE
;  ;          Model_Image_R [x[i],y[i]] = Model_Image_R [x[i],y[i]] + ((loc[4,i])*(atoms_per_packet)*g[i]*1.e-6*(platescale/(body_radius*1.e5))^2.) 
;  ;          ;(fractional packet content)*(atoms per packet)*(photons/(atom*s))*(10^-6 for Rayleighs) * ((pixels/cm)^2)
;  ;        ENDFOR
;  ;        
;    ;optionally draw a point in the sunwardfacing direction at 1 body radii
;    sunward_point = 0
;    IF keyword_set(sunward_point) then begin 
;      cspice_spkpos, 'Sun', ephemeris_time, 'J2000', 'NONE', body, sunward_vector, light_time ;Sun's J2000 position with respect to the body  
;      sunward_vector = sunward_vector / norm(sunward_vector) ;normalize to a unit vector
;      ;IN the "image frame," the +Z axis is the Observer to Body vector and the Sun-Body vector lies in the XZ plane. 
;      cspice_mxv, alignment_matrix, sunward_vector, sunward_vector_prime ;rotate into the "image frame"
;      sun_arrow = platescale * 2. * Body_radius * sunward_vector_prime[0:1] + [x_center, y_center]
;      Model_Image_R[sun_arrow[0],sun_arrow[1]] =  Ceil(Max(Model_Image_R), /L64) 
;    ENDIF

  ;ERROR HANDLING: Make sure no there are no infinite, NaN or negative pixels
  check_bad_pixels = WHERE(FINITE(Model_Image_R , /NAN),count_bad_pixels) 
  if (min(Model_Image_R ) lt 0) OR (fix(n_elements(Model_Image_R )) ne fix(total(finite(Model_Image_R )))) then stop
  save, Model_Image_R , filename = strcompress(directory+'Model_Image_R'+string(loop_number)+'.sav')
  save, Model_Image_CD , filename = strcompress(directory+'Model_Image_CD'+string(loop_number)+'.sav')
endif

if part eq 2 then begin ;when the model is on the last loop plot the averaged output.
  print, 'Generating image. . . '
  loadgood ;load color table 
  body_ang_radius = 4.27 ;HACK HACK! 

  ; Find the body's radius in Km
    cspice_bodvrd, Body,'RADII', 3, Body_radius 
    Body_radius = Body_radius[0]
    x_center=center_in_frame[0]*Output_Size_In_Pixels[0] & y_center=center_in_frame[1]*Output_Size_In_Pixels[1] & z_center=center_in_frame[2]*Output_Size_In_Pixels[2]  
    platescale = float(float(Output_Size_In_Pixels[0]) / (float(plot_range) * body_radius)) ;plate scale of the output image in PIXELS per KM
    stackedModel_Image_R  = MRDFITS(strcompress(directory + Output_Title + '_R.fit'), 0, /SILENT)
    stackedModel_Image_CD  = MRDFITS(strcompress(directory + Output_Title + '_CD.fit'), 0, /SILENT)
    s = size(stackedModel_Image_R)
  
  ; Compute the projected distance anti-sunward for cases where one object viewed from the center of another object
    if (viewpoint ne 'Above Ecliptic' and viewpoint ne 'Moon Spot') then begin
      cspice_spkpos, body, ephemeris_time, 'J2000', 'NONE', '10', body_position, light_time ; bodys J2000 position with respect to the Sun
      cspice_spkpos, viewpoint, ephemeris_time, 'J2000', 'NONE', '10', observer_position, light_time ; Observers J2000 position with respect to the Sun
      observer_to_body = body_position-observer_position ;A vector from the observer location to the body
      S_P_O = cspice_vsep(body_position, observer_to_body) ;the Sun-body-Observer (STO) angle  
      actual_over_apparent = 1./(sin(!pi-S_P_O)) ;the projection effect: down tail length is longer than that projected on the sky 
      radii_per_pixel_crosstail = plot_range/s[1] 
      radii_per_pixel_downtail = plot_range/s[2]        
      pix_per_radii_crosstail = (1./radii_per_pixel_crosstail) ;(rebinned) to 512x512
      pix_per_radii_downtail = (1./radii_per_pixel_downtail) / actual_over_apparent
    endif
  
  ; Draw the outline of the body and a line at the terminator for this phase angle
    if keyword_set(label_phase) then begin
      phase_angle = cspice_phaseq( ephemeris_time, body, 'Sun', viewpoint, 'LT+S' )
      if keyword_set(debug) then print, 'Phase Angle =', strcompress(string(phase_angle*!radeg) + ' Degrees') 
        points = (2. * !PI / 999.0) * FINDGEN(1000) 
      x_circle = x_center + ang_platescale * COS(points) * body_ang_radius
      y_circle = y_center + ang_platescale * SIN(points) * body_ang_radius
      points = (!PI / 999.0) * FINDGEN(1000) + !pi/2.
      x_phase = x_center + (ang_platescale * Body_ang_radius * COS(phase_angle) * COS(points))
      y_phase = y_center + (ang_platescale * Body_ang_radius * SIN(points))
      stackedModel_Image_R [x_circle,y_circle] = Ceil(Max(stackedModel_Image_R), /L64)   ;draw the outline of the body.  
      stackedModel_Image_R [x_phase,y_phase] = Ceil(Max(stackedModel_Image_R), /L64)     ;draw a line at the terminator of the body. 
      stackedModel_Image_CD [x_circle,y_circle] = Ceil(Max(stackedModel_Image_CD), /L64) ;draw the outline of the body.  
      stackedModel_Image_CD [x_phase,y_phase] = Ceil(Max(stackedModel_Image_CD), /L64)   ;draw a line at the terminator of the body.     
    endif

  ; scale axis in RA and Dec
    step = FOV/4. ;this is the x axis step size between ticks in body radii
    xzero = center_in_frame[0]*s[1]

    TICLABELS, Body_dec_WRT_obs*!radeg - FOV/(2.*3600.), N_ticks+1, fov/(60.*(N_ticks)), declabs, delta = 1
    TICLABELS, Body_RA_WRT_obs*!radeg - (FOV*cos(Body_dec_WRT_obs)/(2.*3600.)), N_ticks+1., (fov*cos(Body_dec_WRT_obs)/(60.*N_ticks)) / 15., RAlabs, /RA, delta = 1

    ;window, 0, xs = Output_Size_In_Pixels[0]*3, ys = Output_Size_In_Pixels[1]*3
    axis_format = {XTicks:N_elements(RALABS)-1, XTickname:RALABS, Yticks:N_elements(DecLABS)-1, Ytickname:DecLABS, CHARSIZE:1.2, TICKLEN:-0.01}
    ;cgimage, Model_Image_R, AXKEYWORDS=axis_format, /Axes, /KEEP_ASPECT_RATIO

  !p.charsize = 2. 
  !p.charthick = 2
  cgPS_open, filename = strcompress(directory + Output_title + '_Column_Density.eps'), /ENCAPSULATED, /NOMATCH
   Log_Display = 1
   !P.font=1
   device, SET_FONT = 'Helvetica Bold', /TT_FONT
   ;Warn if optically thick 
      optically_thick = where(stackedModel_Image_CD gt Line_data.unity_optical_depth, count)
      print, count, ' Pixels are optically thick'   
   ;Scale the column density to be plotted over 1-1000 range
      Emission_Scaling = strmid(Image_type, 10) 
      if Emission_Scaling eq 1.e6 then scale = '1.e11'
      if Emission_Scaling eq 1.e5 then scale = '1.e10'
      if Emission_Scaling eq 1.e4 then scale = '1.e9'
      if Emission_Scaling eq 1.e3 then scale = '1.e8'
      if Emission_Scaling eq 1.e2 then scale = '1.e7'
      if Emission_Scaling eq 1.e1 then scale = '1.e6'
      if Emission_Scaling eq 1.e0 then scale = '1.e5'
      Scale_String = strcompress('x 10!U' + strmid(scale, stregex(scale, 'e') + 1) + '!N')               
      Scaled_Model_Image_CD  = stackedModel_Image_CD  / float(scale)   
   cspice_et2utc, ephemeris_time, "C", 0, utc_current       
      MinValue = Floor(Min(Scaled_Model_Image_CD))
      MaxValue = 1.e3 ;Ceil(Max(Scaled_Model_Image_CD))
      ;format x axis       
        step = plot_range/4. ;this is the x axis step size between ticks in body radii
        xzero = center_in_frame[0]*s[1]
        ; Normal Cases: One object viewed from the center of another object
        if (viewpoint ne 'Above Ecliptic' and viewpoint ne 'Moon Spot') then $
          xtickv = ((step*(findgen(7) - 3))*pix_per_radii_downtail)+xzero $        
          else xtickv = step*(findgen(7) - 3)+xzero 
        xticks = n_elements(xtickv)-1
        xtickname = step*(findgen(7) - 3)
        xtickname = strcompress(string(xtickname, format = '(F10.2)'),/remove_all)       
      ;format y axis  
        step = plot_range/5. ;this is the y axis step size between ticks in body radii
        yzero = center_in_frame[1]*s(2)
        if (viewpoint ne 'Above Ecliptic' and viewpoint ne 'Moon Spot') then $
          ytickv = step*(findgen(7) - 3)*pix_per_radii_crosstail+yzero $
          else ytickv = step*(findgen(7) - 3)+yzero 
        yticks = n_elements(ytickv)-1
        ytickname = step*(findgen(7) - 3)
        ytickname = strcompress(string(ytickname, format = '(F10.2)'),/remove_all) 
      if not keyword_set(Log_Display) then begin
        loadct, 0, /silent
        cgImage, Scaled_Model_Image_CD, Stretch=1, MinValue = MinValue, $
          MaxValue = MaxValue, Position = [0.125, 0.125, 0.9, 0.800], /keep_aspect, $
          oposition = test, Title = strcompress('Model result - ' + strmid(utc_current,0,17)), charsize =1.2, /Window, /save
        AXIS,0,0,0, YAXIS = 0, YSTYLE = 1, color=0., charsize=1.3,ythick=3.,$
            yticks = yticks, ytickv = ytickv, ytickname = ytickname, yticklen = -.01, $
            Ytitle = STRCOMPRESS(Body + ' Radii'), /data 
        AXIS,0,0,0, XAXIS = 0, XSTYLE = 1, color=0., charsize=1.3,Xthick=3.,$
            xticks = xticks, xtickv = xtickv, xtickname = transpose(xtickname), xticklen = -.01, /data, $
            Xtitle = strcompress('Sunward Distance ' + body + ' Radii')
        cgColorbar, /Vertical, Position = [!x.window[1]+.07, !y.window[0], !x.window[1]+.1, !y.window[1]], Range=[MinValue, MaxValue], $
          Title = strcompress('!3' + Particle_data.name + ' Column Density ('+ scale_string +' cm!U-2!N)'), TLocation='Left' ,/right             
      endif else begin   ; Log scaled image 
        loadct, 0, /SILENT
        number_of_colors = 8                                                    ;number of colors / labels on the color bar, loadgood has 8
        colorbar_label = reverse(MaxValue*.5^findgen(number_of_colors + 1.))    ;color bar labels
        colorbar_label[0] = 0                                                   ;force the bottom label to zero, true zero has a negaitve inf. scaled value 
        colorbar_label = strtrim(string(colorbar_label, format = '(i)'), 2)     ;display color bar labels as integers
        ;Scale the data so that factors of 2 represents integer increases in brightness
        Log_scale_img_CD = alog(float(Scaled_Model_Image_CD) / MaxValue) / alog(2.) + number_of_colors
        ;number_of_colors is the value of the new image at the MaxValue peak level
        cgImage, Log_scale_img_CD, Stretch=1, MinValue = MinValue, MaxValue = number_of_colors, Position = [0.125, 0.125, 0.9, 0.800], /keep_aspect, $
          oposition = test, Title = strcompress('Model result - ' + strmid(utc_current,0,17)), charsize =1.2, /Window, /save
        !p.charsize = 2. 
        AXIS,0,0,0, YAXIS = 0, YSTYLE = 1, color=0., charsize=1.3, ythick=3.,$
            yticks = yticks, ytickv = ytickv, ytickname = ytickname, yticklen = -.01, $
            Ytitle = STRCOMPRESS(Body + ' Radii'), /data 
        AXIS,0,0,0, XAXIS = 0, XSTYLE = 1, color=0., charsize=1.3, Xthick=3.,$
            xticks = xticks, xtickv = xtickv, xtickname = transpose(xtickname), xticklen = -.01, /data, $
            Xtitle = strcompress('Sunward Distance ' + body + ' Radii')
        cgColorbar, /Vertical, Position = [!x.window[1]+.07, !y.window[0], !x.window[1]+.1, !y.window[1]], Range=[MinValue, MaxValue], $
          Title = strcompress('!3' + Particle_data.name + ' Column Density ('+ scale_string +' cm!U-2!N)'), TLocation='Left', /right, $  
          ticknames = colorbar_label, divisions = number_of_colors, ticklen = 0, charsize = 1.5
      
        if keyword_set(Label_Time) then begin   
            xyouts, .66,.75, strcompress(Label_Time + ' Hours'), /normal, color=cgColor('white'), charsize = 1.5           
        endif     
      endelse           
  cgPS_Close 

  cgPS_Open, filename=strcompress(directory + Output_title + '_Emission.eps'), /ENCAPSULATED, /NOMATCH
  Log_Display = 1 
  !P.font=1
  device, SET_FONT = 'Helvetica Bold', /TT_FONT
      Emission_Scaling = strmid(Image_type, 10) 
      if Emission_Scaling eq 1.e6 then Emission_Label = 'MR'
      if Emission_Scaling eq 1.e5 then Emission_Label = 'x100 kR'
      if Emission_Scaling eq 1.e4 then Emission_Label = 'x10 kR'
      if Emission_Scaling eq 1.e3 then Emission_Label = 'kR'
      if Emission_Scaling eq 1.e2 then Emission_Label = 'x100 R'
      if Emission_Scaling eq 1.e1 then Emission_Label = 'x10 R'
      if Emission_Scaling eq 1.e0 then Emission_Label = 'R'        
      Scaled_Model_Image_R  = stackedModel_Image_R / Emission_Scaling ;convert to the Model_Image_R  to match the color bar scale  
      if (viewpoint eq 'STEREO AHEAD') or (viewpoint eq 'STEREO BEHIND') then begin
        loadct,0, /SILENT
        convolve_kernal = psf_Gaussian( NPIXEL = 10., FWHM = [1.5,1.5], /NORMALIZE) ;see Bewsher et al. 2010 THESE ARE APPROXIMATE PSFs!
        Scaled_Model_Image_R  = convol(Scaled_Model_Image_R ,convolve_kernal, /normalize);convolve it with a stereo fit PSF
        Scaled_Model_Image_R  = rebin(Scaled_Model_Image_R ,512,512) ;rebin it to match the stereo size
      endif     
      !p.charsize=1.5
      cspice_et2utc, ephemeris_time, "C", 0, utc_current       
      minValue = Floor(Min(Scaled_Model_Image_R))
      maxValue = 1.e3 ;Ceil(Max(Scaled_Model_Image_R))
      ;format x axis       
        step = plot_range/4. ;this is the x axis step size between ticks in bodyary radii
        xzero = center_in_frame[0]*s(1)
        if (viewpoint ne 'Above Ecliptic' and viewpoint ne 'Moon Spot') then $
          xtickv = step*(findgen(7) - 3)*pix_per_radii_downtail+xzero $        
          else xtickv = step*(findgen(7) - 3)+xzero          
        xticks = n_elements(xtickv)-1
        xtickname = step*(findgen(7) - 3)
        xtickname = strcompress(string(xtickname, format = '(F10.2)'),/remove_all)         
      ;format y axis  
        step = plot_range/5. ;this is the y axis step size between ticks in bodyary radii
        yzero = center_in_frame[1]*s(2)
        if (viewpoint ne 'Above Ecliptic' and viewpoint ne 'Moon Spot') then $
          ytickv = step*(findgen(7) - 3)*pix_per_radii_crosstail+yzero $
          else ytickv = step*(findgen(7) - 3)+yzero 
        yticks = n_elements(ytickv)-1
        ytickname = [step*(-3.),step*(-2.),step*(-1.),0.,step*1.,step*2.,step*3.]
        ytickname = strcompress(string(ytickname, format = '(F10.2)'),/remove_all)                   
      if not keyword_set(Log_Display) then begin
        loadct,3, /SILENT
        cgImage, Scaled_Model_Image_R, Stretch=1, MinValue = Floor(Min(Scaled_Model_Image_R)), $
          MaxValue = Ceil(Max(Scaled_Model_Image_R)), Position = [0.125, 0.125, 0.9, 0.800], /keep_aspect, $
          oposition = test, Title = strcompress('Model result - ' + strmid(utc_current,0,17)), charsize =1.2, /Window, /save
        AXIS,0,0,0, YAXIS = 0, YSTYLE = 1, color=0., charsize=1.3,ythick=3.,$
          yticks = yticks, ytickv = ytickv, ytickname = ytickname, yticklen = -.01, $
          Ytitle = STRCOMPRESS(Body + ' Radii'), /data 
        AXIS,0,0,0, XAXIS = 0, XSTYLE = 1, color=0., charsize=1.3,Xthick=3.,$
          xticks = xticks, xtickv = xtickv, xtickname = transpose(xtickname), xticklen = -.01, /data, $
          Xtitle = strcompress('Sunward Distance ' + body + ' Radii')
        cgColorbar, /Vertical, Position = [!x.window[1]+.07, !y.window[0], !x.window[1]+.1, !y.window[1]], Range=[MinValue, MaxValue], $
          Title = Strcompress('Simulated ' + Line_data.name + ' Line Emission ('+Emission_Label+')'), TLocation='Left', /right, charsize = 1.5              
      endif else begin   ; Log scaled image
        loadct,3, /SILENT
        ;loadgood  
        number_of_colors = 8                                                    ;number of colors / labels on the color bar, loadgood has 8
        colorbar_label = reverse(MaxValue*.5^findgen(number_of_colors + 1.))    ;color bar labels
        colorbar_label[0] = 0                                                   ;force the bottom label to zero, true zero has a negaitve inf. scaled value 
        colorbar_label = strtrim(string(colorbar_label, format = '(i)'), 2)     ;display color bar labels as integers
        ;Scale the data so that factors of 2 represents integer increases in brightness
        Log_scale_img_R = alog(float(Scaled_Model_Image_R) / MaxValue) / alog(2.) + number_of_colors
        ;number_of_colors is the value of the new image at the MaxValue peak level
        cgImage, Log_scale_img_R, Stretch=1, MinValue=MinValue, MaxValue = number_of_colors, Position = [0.125, 0.125, 0.9, 0.800], /keep_aspect, $
          oposition = test, Title = strcompress('Model result - ' + strmid(utc_current,0,17)), charsize =1.2, /Window, /save
        AXIS,0,0,0, YAXIS = 0, YSTYLE = 1, color=0., charsize=1.3,ythick=3.,$
          yticks = yticks, ytickv = ytickv, ytickname = ytickname, yticklen = -.01, $
          Ytitle = STRCOMPRESS(Body + ' Radii'), /data 
        AXIS,0,0,0, XAXIS = 0, XSTYLE = 1, color=0., charsize=1.3,Xthick=3.,$
          xticks = xticks, xtickv = xtickv, xtickname = transpose(xtickname), xticklen = -.01, /data, $
          Xtitle = strcompress('Sunward Distance ' + body + ' Radii')
        cgColorbar, /Vertical, Position = [!x.window[1]+.07, !y.window[0], !x.window[1]+.1, !y.window[1]], Range=[MinValue, MaxValue], $
          Title = Strcompress('Simulated ' + Line_data.name + ' Line Emission ('+Emission_Label+')'), TLocation='Left', /right, $
          ticknames = colorbar_label, divisions = number_of_colors, ticklen = 0, charsize = 1.5 
      endelse

      if keyword_set(Label_North) then begin   
        ;recall the Northward unit vector and transform it from the J2000 frame to the "tail aligned frame"
        restore, strcompress(directory+'North_tail_aligned_frame.sav'), DESCRIPTION = North_vector_info 
         if strcompress(STRMID(North_vector_info, 46),/remove_all) eq strcompress(ephemeris_time,/remove_all) then begin              
            if North_tail_aligned_frame[1] gt 0 then arrow, .715, .82, .715 + North_tail_aligned_frame[0]/12., $
              .82 + North_tail_aligned_frame[1]/12., /normalized, color=cgColor('dodger blue')   , thick = 5., hthick = 5., HSIZE = 500 $ 
            else arrow, .715, .9, .715 + North_tail_aligned_frame[0]/12., $
              .90 + North_tail_aligned_frame[1]/12., /normalized, color=cgColor('dodger blue')   , thick = 5., hthick = 5., HSIZE = 500 
            xyouts, .6,.85, 'NORTH',/normal, color=cgColor('dodger blue')          
         endif else begin 
            print, 'Error: Time of northward vector definition does not match time simulated'   
            cspice_et2utc, double(strcompress(STRMID(North_vector_info, 46))), "C", 0, utc_North_Defined
            print, 'North was defined for ', utc_North_Defined
            print, 'Time being simulated is ', utc_current
            stop 
         endelse
      endif  
      
      if keyword_set(Label_Time) then begin   
            xyouts, .66,.75, strcompress(Label_Time + ' Hours'), /normal, color=cgColor('white'), charsize = 1.5          
      endif             
      cgPS_Close 
endif
;*************************************************************************************************************************************************************
return

end