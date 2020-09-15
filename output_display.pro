pro output_display, loc, g, atoms_per_packet, Image_type, loop_number, part, Label_North = Label_North, Label_phase = Label_phase, $
    Label_time = Label_time 

; REVISION HISTORY:
;    Written:                     C. Schmidt, July 2009
;    -Vectorized some, added 'Above Ecliptic' & 'Moon Spot' viewpoints
;    -Changed Plot_Range to FOV: from planetary radii to arcseconds 
;                                 C. Schmidt, Nov 2018

COMMON Output_shared, Plot_range, Output_Size_In_Pixels, Output_Title, Center_in_frame, viewpoint, FOV, N_ticks, Tickstep, Observatory, Above_Ecliptic, Boresight_Pixel, Aperture_Corners
COMMON Model_shared, Body, Ephemeris_time, Seed, Directory, Particle_data, Line_data, Debug

;stretch = stretch                  ; Display image in log scaling
loadct, 3, /SILENT                  ; load color table, use color table 3 for a sodium look, or loadgood.pro for the BU stepped rainbow color table

N_particles = N_elements(loc[0,*])  ; Number of particles simulated

; ------------------------------------------------- Determine the RA, DEC, & ANGULAR SIZE of the Body's disk -------------------------------------------------------

; Find the body's radius in Km
  cspice_bodvrd, Body, 'RADII', 3, Body_radius
  Body_radius = Body_radius[0]

  Case viewpoint of
    'Earth' or 'Moon Spot': begin

      ; Get the body-fixed non-inertial IAU_Earth Coordinates of the Observatory in Cartesian
        cspice_bodvrd, 'EARTH', 'RADII', 3, radii
        flat = (radii[0] - radii[2])/radii[0]
        OBSERVATORY_cs, observatory, obs_struct ;longitude is in degrees *West*, need to switch to East Longit for J2000 conversion via cspice_georec
        cspice_georec, (360. - obs_struct.longitude) * CSPICE_RPD(), obs_struct.latitude * CSPICE_RPD(), obs_struct.altitude / 1.e3, radii[0], flat, obs_IAU_Earth
          
      ; Convert to J2000
        cspice_pxform, 'IAU_EARTH', 'J2000', ephemeris_time, Earth_to_J2000_transform_matrix 
        obs_J2000 = TRANSPOSE(Earth_to_J2000_transform_matrix) # obs_IAU_Earth
    
      ; Get the RA and Dec of the Body 
        cspice_spkpos, Body, ephemeris_time, 'J2000', 'LT', 'Earth', body_Earth_vector, light_time 
        cspice_recrad, (body_Earth_vector - obs_J2000), Body_distance_WRT_obs, Body_RA_WRT_obs, Body_dec_WRT_obs
        ; VALIDATED Above result: Moon's RA and Dec from an observatory accurate within ~3" via JPL Horizons.
    
      ; Get the angular radius of the Body in arcseconds          
        body_ang_radius = asin(Body_radius / Body_distance_WRT_obs)*3600.d/CSPICE_RPD()
      end
    else: begin
     
      ; Get the RA and Dec of the Body
        cspice_spkpos, Body, ephemeris_time, 'J2000', 'LT', Viewpoint, body_Obs_vector, light_time
        cspice_recrad, body_obs_vector, Body_distance_WRT_obs, Body_RA_WRT_obs, Body_dec_WRT_obs
  
      ; Get the angular radius of the Body in arcseconds
        body_ang_radius = asin(Body_radius / Body_distance_WRT_obs)*3600.d/CSPICE_RPD()
      end
  endcase

;------------------------------------------- For the Moon Spot, we'll also want the ANTI-SOLAR RA & DEC -----------------------------------------------------------
  Case viewpoint of 
    'Moon Spot': begin
    ; Get the RA and Dec of the anti-sunward vector Sun WRT the observer
      cspice_spkpos, 'Sun', ephemeris_time, 'J2000', 'LT', 'Earth', Sun_Earth_vector, light_time 
      cspice_recrad, -(Sun_Earth_vector-obs_J2000), radius, anti_Sun_RA, anti_Sun_dec
      cspice_recrad, (Sun_Earth_vector-obs_J2000), radius, Sun_RA, Sun_dec

      cspice_spkpos, Body, ephemeris_time, 'J2000', 'LT', 'Earth', Body_Earth_vector, light_time 
      cspice_recrad, -(Body_Earth_vector-obs_J2000), radius, anti_Body_RA, anti_Body_dec
      cspice_recrad, (Body_Earth_vector-obs_J2000), radius, Body_RA, Body_dec

      cspice_spkpos, 'Jupiter', ephemeris_time, 'J2000', 'LT', 'Earth', Jupiter_Earth_vector, light_time 
      cspice_recrad, (Jupiter_Earth_vector-obs_J2000), radius, Jup_RA, Jup_dec
      
    ; Now get the RA and Dec of all the particles at the "image taken" ephemeris time, when the integration is done
      ; Test particles are in the J2000 frame with body/Lunar center as the origin 
      ; Move Loc array origin to observatory, still in J2000
      cspice_spkpos, body, ephemeris_time, 'J2000', 'LT', 'Earth', Body_Earth_vector, light_time        
      Loc_WRT_obs_J2000 = loc[0:2,*] + rebin(Body_Earth_vector, 3, N_particles) - rebin(obs_J2000, 3, N_particles) 
      cspice_recrad, Loc_WRT_obs_J2000, Loc_radius_WRT_obs, Loc_RA_WRT_obs, Loc_dec_WRT_obs

;        ;looking sunward
;        window, 0, xs = 500, ys = 500
;;        x = (Loc_RA_WRT_obs - Sun_RA)/CSPICE_RPD()
;;        Y = (Loc_dec_WRT_obs - Sun_dec)/CSPICE_RPD()
;        x = (Loc_RA_WRT_obs - Body_RA)/CSPICE_RPD()
;        Y = (Loc_dec_WRT_obs - Body_dec)/CSPICE_RPD()
;        cgplot, x, y, xr = [-1,1], yr = [-1,1], psym = 4    
;        
;        ;looking anti-sunward
;        window, 1, xs = 500, ys = 500
;        x = (Loc_RA_WRT_obs - anti_Sun_RA)*3600.d/CSPICE_RPD()   ; Units of Arcseconds
;        Y = (Loc_dec_WRT_obs - anti_Sun_dec)*3600.d/CSPICE_RPD() ; Units of Arcseconds
;        cgplot, x, y, xr = [-20, 20]*3600.d, yr = [-20, 20]*3600.d, psym = 4    
;        stop    
      end
    else:        
  endcase     

;---------------------------FROM THE OBSERVER'S VIEWPOINT: Determine which particles are occulted by the body------------------------

  ; Find the rays that intersect the body  
  ; Then indentify the subset of these particles that are behind the planet
  ; Particles behind the planet cannot emit. Particles in front can emit.
    occulted = replicate(1., N_particles)  ; set occulted elements to 0. This array will be multiplied in this program's column density and brightness calculation  
    cspice_spkpos, body, ephemeris_time, 'J2000', 'LT', viewpoint, Body_observer_vector, light_time
    J2000_rays = loc[0:2,*] + rebin(Body_observer_vector, 3, N_particles)
    for i  = 0, N_particles-1 do begin
      if loc[4,i] eq 0. then continue ; these particles are not in the exosphere and have already settled onto the surface
      
      cspice_sincpt, 'Ellipsoid', body, Ephemeris_time, 'IAU_'+body, 'LT', $
                      Viewpoint, 'J2000', J2000_rays[*,i], spoint, trgepc, srfvec, found
      ; if no surface intercept is found, the ray to this particle does not intersect the planet
      if (found eq 1) and (norm(J2000_rays[*,i]) gt norm(srfvec)) then occulted[i] = 0. ; this particle is behind the target body
    endfor
  
  Print, 'This section should be tested against various particle release lat/lon locations' 

;----------------------This section builds a 2D image from a specified look direction------------------

; This performed in an angular sense, pixels beween xi & xj and yi & yj angles from a center axis fall into pixel xi and xj 
; An exception is the 'Above Ecliptic' viewpoint perspective, which creates a loc_prime array where the z axis is normal to the ecliptic in positive north. 

ang_platescale = float(Output_Size_In_Pixels[0]) / float(FOV) ;plate scale of the output image in PIXELS PER ARCSEC 

if part eq 1 then begin
  
  ; Looking down over the ecliptic is a special case
    if Keyword_set(Above_ecliptic) then begin  
        cspice_pxform, "J2000", "GSE", ephemeris_time, J2000_to_GSE_transform_matrix ;'Z' axis is normal to the ecliptic now.
        loc_prime_above_ecliptic = TRANSPOSE(J2000_to_GSE_transform_matrix) # loc[0:2,*]        
        radii_per_pixel = plot_range / Output_Size_In_Pixels[0] ;HACK HACK HACK NEED X AND Y IN CASE IT'S NOT CENTERED IN THE FRAME
    endif
    
  ; Normal Cases: One object viewed from the center of another object
    if (viewpoint ne 'Above Ecliptic' and viewpoint ne 'Moon Spot') then begin
      
      ; ********************************************Loc_Prime Frame Definition***********************************************************
      ; *                                                                                                                               *
      ; *  Define a coordinate system in which the +Z axis is the Observer to Body vector and the Sun-Body vector lies in the XZ plane. *
      ; *  This means that in Loc-prime, the projection of the sunward direction is always in the -X dimension                          *
      ; *                                                                                                                               *
      ; *********************************************************************************************************************************
      
      ; Generate a coordinate system where the +Z axis lies along the Observer-Body line, and the *projected* anti-sunward tail is along the +X axis.
        cspice_spkpos, 'Sun', ephemeris_time, 'J2000', 'NONE', Body, Body2Sun_vector, light_time
        cspice_spkpos, body, ephemeris_time, 'J2000', 'NONE', viewpoint, Observer2Body_vector, light_time
    
      ; Print the phase angle and distance 
        if Keyword_Set(debug) then begin
          print, body, ' phase angle from ', viewpoint, ':',  !radeg*cspice_vsep( Observer2Body_vector, Body2Sun_vector )
          print, viewpoint, ' distance to ', body, ' center:', norm(Observer2Body_vector), ' km'
        endif
          
      ; Now that we have the Observer to Body and Sun to Body vectors, we can define a new frame for visualization
      ; Find the rotation matrix to transform particle locations from their current J2000 Body-centered frame into the desired loc_prime coordinate frame.
        cspice_twovec, Observer2Body_vector, 3, Body2Sun_vector, 1, alignment_matrix  ;Re-defines Observer2Body_vector as the '+Z' axis
      
      ; Put the particle positions into this new "aligned" coordinate frame: loc_prime
        loc_prime = TRANSPOSE(alignment_matrix) # loc[0:2,*]
    endif  

  ;------------------------------------------------BUILD AN IMAGE VIEWED FROM A FINITE RADIUS USING AN ANGULAR PIXEL SIZE-----------------------------------------------------------------
     ; Define the empty image arrays that this loop will write into 
       Model_Image_R  = fltarr(Output_Size_In_Pixels[0], Output_Size_In_Pixels[1]) ;Image in Rayleighs in the XY plane to be written 
       Model_Image_CD = fltarr(Output_Size_In_Pixels[0], Output_Size_In_Pixels[1]) ;Image in cgs column density in the XY plane to be written 
     
     ; Place the body center in pixel space
       x_center = center_in_frame[0]*Output_Size_In_Pixels[0]
       y_center = center_in_frame[1]*Output_Size_In_Pixels[1]
       z_center = center_in_frame[2]*Output_Size_In_Pixels[2]

   ; Get the RA and Dec of the 'Body' and all the particles from the observers location...  
     Case viewpoint of
       'Moon Spot' or 'Earth': begin  ; Earth based observations should include the parallax of the topocentric obsever location
          
        ; Get the body-fixed non-inertial IAU_Earth Coordinates of the Observatory in Cartesian
          cspice_bodvrd, 'EARTH', 'RADII', 3, radii
          flat = (radii[0] - radii[2])/radii[0]
          OBSERVATORY_cs, observatory, obs_struct ;longitude is in degrees *West*, need to switch to East Longit for J2000 conversion via cspice_georec
          cspice_georec, (360. - obs_struct.longitude) * CSPICE_RPD(), obs_struct.latitude * CSPICE_RPD(), obs_struct.altitude / 1.e3, radii[0], flat, obs_IAU_Earth
      
        ; convert to J2000
          cspice_pxform, 'IAU_EARTH', 'J2000', ephemeris_time, Earth_to_J2000_transform_matrix 
          obs_J2000 = TRANSPOSE(Earth_to_J2000_transform_matrix) # obs_IAU_Earth

        ; Get the RA and Dec of the Body
          cspice_spkpos, Body, ephemeris_time, 'J2000', 'LT', viewpoint, Body_Observer_vector, light_time
          cspice_recrad, (Body_Observer_vector - obs_J2000), Body_distance_WRT_obs, Body_RA_WRT_obs, Body_dec_WRT_obs
          ;RA and Dec from an observatory verified accurate within ~0.055" via JPL Horizons (for Sun & Mercury)

        ; Now get the RA and Dec of all the particles
          ; All test particles are currently in the J2000 frame with body center as the origin. 
          ; Move origin to observer and convert to lat, lon (RA & Dec)
          cspice_recrad, loc[0:2,*] + rebin(Body_Observer_vector, 3, N_particles) - rebin(obs_J2000, 3, N_particles), $
                         Loc_distance_WRT_obs, Loc_RA_WRT_obs, Loc_dec_WRT_obs
        
        ; Find the Delta RA and Delta Dec FROM BODY CENTER
          if viewpoint eq 'Earth' then begin
            delta_RA        = (Loc_RA_WRT_obs - Body_RA_WRT_obs)*3600.d/CSPICE_RPD() ; Units of Arcseconds
            delta_Dec       = (Loc_dec_WRT_obs-Body_dec_WRT_obs)*3600.d/CSPICE_RPD() ; Units of Arcseconds
          endif
       
        ; For viewing the Moon Spot, this Delta_RA and Delta_Dec will instead defined FROM THE ANTI-SOLAR POINT
          if viewpoint eq 'Moon Spot' then begin
            ;delta_RA  = (Loc_RA_WRT_obs - Jup_RA)*3600.d/CSPICE_RPD() ; Units of Arcseconds
            ;delta_Dec = (Loc_dec_WRT_obs - Jup_dec)*3600.d/CSPICE_RPD() ; Units of Arcseconds
            delta_RA  = (Loc_RA_WRT_obs - anti_Sun_RA)*3600.d/CSPICE_RPD() ; Units of Arcseconds
            delta_Dec = (Loc_dec_WRT_obs-anti_Sun_dec)*3600.d/CSPICE_RPD() ; Units of Arcseconds
          endif             
          window, 0, xs = 500, ys = 500
          plot, delta_RA, delta_Dec, psym = 3
      end
    else: begin ; RA and Dec coordinates make little sense for non-Earth-based observers
      ; Get the RA and Dec of the Body
        cspice_spkpos, Body, ephemeris_time, 'J2000', 'LT', viewpoint, Body_Observer_vector, light_time
        cspice_recrad, Body_Observer_vector, Body_distance_WRT_obs, Body_RA_WRT_obs, Body_dec_WRT_obs
  
      ; Now get the ANGULAR distance between all particles
        cspice_recrad, loc_prime[0:2,*] + rebin([0,0,Body_distance_WRT_obs], 3, N_particles), $
                       Loc_distance_WRT_obs, Loc_RA_WRT_obs, Loc_dec_WRT_obs
      
        angular_radius = abs(Loc_dec_WRT_obs - !pi/2.)
        delta_RA  = angular_radius * cos(Loc_RA_WRT_obs) * 3600.d/CSPICE_RPD() ;this is just an X, Y coordinate in arcseconds from body center
        delta_Dec = angular_radius * sin(Loc_RA_WRT_obs) * 3600.d/CSPICE_RPD()
      end
    endcase    

    ;      window, 0, xs = 500, ys = 500
    ;      plot, delta_RA, delta_Dec, psym = 3
    ;      window, 1, xs = 500, ys = 500
    ;      plot, loc_prime[0,*], loc_prime[1,*], psym = 3

    ; Find the x & y pixel coordinates of the Body's disk
      points   = FINDGEN(1000) * 2. * !PI / 999. 
      x_circle = x_center + ang_platescale * COS(points) * body_ang_radius
      y_circle = y_center + ang_platescale * SIN(points) * body_ang_radius
      disc = polyfillv(x_circle, y_circle, Output_Size_In_Pixels[0], Output_Size_In_Pixels[1])
      if disc[0] ne -1 then begin ;if the plate scale is such that a pixel fits on the bodys disc, get the coordinates of the disc.   
        body_resolved = 'true'                                    ; Signifies there will be pixels on the body's disk 
        dims = SIZE(Model_Image_R, /DIMENSIONS) 
        disc_coordinates = array_indices(dims, disc, /DIMENSIONS) ; Coordinates of pixels occuped by body's disc in x and y
      endif else body_resolved = 'false'    
           
    ; Build up the 2D image. Put the packets into the pixel bins. 
      for i = 0, N_particles-1 do begin 
        if loc[4,i] eq 0. then continue              ; Skip particles that have been fully lost via surface settling
        x = -delta_RA[i] * ang_platescale + x_center ; RA increases **EASTWARD** not westward in the sky
        y = delta_Dec[i] * ang_platescale + y_center 

;        ; Check if the packet's location is on the disc. If so, only display packets between the body and the observer 
;          ;if body_resolved then here = where(disc_coordinates[0,*] eq fix(x) and disc_coordinates[1,*] eq fix(y), count) else count = 0.
;          if body_resolved then here = where(disc_coordinates[0,*] eq round(x) and disc_coordinates[1,*] eq round(y), count) else count = 0.
;          IF count NE 0 THEN on_disc=1 ELSE on_disc=0 
;          if (on_disc eq 1) and (Loc_distance_WRT_obs[i] gt Body_distance_WRT_obs) then continue ;skip packets behind the disc      
    
        ; Allocate x,y positions to pixels in the 2d Model_Image_R if within the plot window 
          if x ge 0 and x lt Output_Size_In_Pixels[0] and y ge 0 and y lt Output_Size_In_Pixels[1] then begin 
            
            ; Output an image of the column density
              Model_Image_CD[x,y] = Model_Image_CD[x,y] + (loc[4,i] * occulted[i] * atoms_per_packet * ( ang_platescale*body_ang_radius/(1.e5*body_radius) )^2.)
              ; atoms / cm^2 = (fractional packet content) * (particle is visible?) * (atoms per packet) * ((pixels/cm)^2)  
              ; where pixels/cm = (pixels/arcsec) * (arcsec/cm)
          
            ; Output an image in rayleighs
              Model_Image_R[x,y] = Model_Image_R[x,y] + (loc[4,i] * occulted[i] * atoms_per_packet * g[i] * 1.e-6 * ( ang_platescale*body_ang_radius/(1.e5*body_radius) )^2.)              
              ;(fractional packet content) * (particle is visible?) * (atoms per packet) * (photons/(atom*s)) * (10^-6 for Rayleighs) * ((pixels/cm)^2)   
          endif
      endfor 
      
      ;window, 0
      ;cgimage, Model_Image_CD, /keep_aspect_ratio
;stop
  ;------------------------------------------BUILD AN IMAGE VIEWED FROM INFINITY LOOKING DOWN FROM ABOVE THE ECLIPTIC--------------------------------------------------------
if Keyword_set(Above_ecliptic) then begin  
    x_center=center_in_frame[0]*Output_Size_In_Pixels[0] & y_center=center_in_frame[1]*Output_Size_In_Pixels[1] & z_center=center_in_frame[2]*Output_Size_In_Pixels[2]  
    platescale = float(float(Output_Size_In_Pixels[0]) / (float(plot_range) * body_radius)) ;plate scale of the output image in PIXELS per KM
    Above_Ecliptic_Image_R  = fltarr(Output_Size_In_Pixels[0], Output_Size_In_Pixels[1]) ;Image in Rayleighs in the XY plane to be written 
    Above_Ecliptic_Image_CD = fltarr(Output_Size_In_Pixels[0], Output_Size_In_Pixels[1]) ;Image in cgs column density in the XY plane to be written 
    ;Above_Ecliptic_Image_V  = fltarr(Output_Size_In_Pixels[0], Output_Size_In_Pixels[1]) ;Image in km/s to be written
    ;optically_thick = Above_Ecliptic_Image_R  ;use this to generate a seperate image for optically thick pixels, then replace where needed in the image.      
    
  ; Get the coordinates in the x-y plane for the body's disc
    points   = (2. * !PI) * FINDGEN(1000) / 999.0
    x_circle = x_center + platescale * COS(points) * body_radius
    y_circle = y_center + platescale * SIN(points) * body_radius
    disc = polyfillv(x_circle, y_circle, Output_Size_In_Pixels[0], Output_Size_In_Pixels[1])
    if disc[0] ne -1 then begin ;if the plate scale is such that a pixel fits on the bodys disc, get the coordinates of the disc.   
      body_resolved = 1. 
      dims = SIZE(Above_Ecliptic_Image_R, /DIMENSIONS) 
      disc_coordinates = array_indices(dims, disc, /DIMENSIONS) ;coordinates of pixels occupes by Body's disc in x and y
    endif else body_resolved = 0.     
    
    ; Put the packets into the pixel bins
      for i = 0, N_particles-1 do begin 
        if (loc[4,i] eq 0.) then continue ;To save computation time, skip particles that have been lost via surface collisions
        
        x = loc_prime_above_ecliptic[0,i] * platescale + x_center ;(distance from Body in km) * (pixels per km) + (Body center pixel location, x) 
        y = loc_prime_above_ecliptic[1,i] * platescale + y_center ;(distance from Body in km) * (pixels per km) + (Body center pixel location, y) 
        z = loc_prime_above_ecliptic[2,i] * platescale + z_center ;(distance from Body in km) * (pixels per km) + (Body center pixel location, z)      
        ;      v_x = loc_prime[5,i] ; km/s
        ;      v_y = loc_prime[6,i] ; km/s
        ;      v_z = loc_prime[7,i] ; km/s
    
        ;if the packet location is on the disc show only packets in the negative z (between the body and the observer) 
        if body_resolved then here = where(disc_coordinates[0,*] eq fix(x) and disc_coordinates[1,*] eq fix(y), count) else count = 0.
        IF count NE 0 THEN on_disc=1. ELSE on_disc=0. 
        if (on_disc eq 1) and (loc_prime_above_ecliptic[2,i] lt 0.) then continue ;skip packets behind the disc      
          if x ge 0 and x lt Output_Size_In_Pixels[0] and y ge 0 and y lt Output_Size_In_Pixels[1] then begin ;allocate x,y positions to pixels in the 2d Above_Ecliptic_Image_R  if within the plot window 
              ;output an image of the column density
              Above_Ecliptic_Image_CD[x,y] = Above_Ecliptic_Image_CD[x,y] + ((loc[4,i]) * (atoms_per_packet) * (platescale/1.e5)^2.)
              ;(fractional packet content) * (atoms per packet) * ((pixels/cm)^2)  
            
              ;output an image in rayleighs
              Above_Ecliptic_Image_R[x,y] = Above_Ecliptic_Image_R[x,y] + ((loc[4,i]) * (atoms_per_packet) * g[i] * 1.e-6 * (platescale/1.e5)^2.)              
              ;(fractional packet content) * (atoms per packet) * (photons/(atom*s)) * (10^-6 for Rayleighs) * ((pixels/cm)^2)            
           endif
      endfor
endif ;----when next editing: try to combine this with the above regular loop!

  ; ERROR HANDLING: Make sure no there are no infinite, NaN or negative pixels
    check_bad_pixels = WHERE( FINITE(Model_Image_R, /NAN), count_bad_pixels) 
    if count_bad_pixels gt 0 then stop
    if (min(Model_Image_R ) lt 0) OR (fix(n_elements(Model_Image_R)) ne fix(total(finite(Model_Image_R)))) then stop
  
  ; Write the output frames to save files
    save, Model_Image_R , filename = strcompress(directory+'Model_Image_R'+string(loop_number)+'.sav')
    save, Model_Image_CD, filename = strcompress(directory+'Model_Image_CD'+string(loop_number)+'.sav')
    if Keyword_set(Above_ecliptic) then save, Above_Ecliptic_Image_R , filename = strcompress(directory+'Above_Ecliptic_Image_R'+string(loop_number)+'.sav')
    if Keyword_set(Above_ecliptic) then save, Above_Ecliptic_Image_CD, filename = strcompress(directory+'Above_Ecliptic_Image_CD'+string(loop_number)+'.sav')
endif

if part eq 2 then begin ;Once the model reaches it's last loop, plot the averaged output.
  print, 'Displaying image. . . '
  ; Restore the images that where written in the main routine as an average over all the loops 
    stackedModel_Image_R  = MRDFITS(strcompress(directory + Output_Title + '.fit'), 0, /SILENT)
    stackedModel_Image_CD = MRDFITS(strcompress(directory + Output_Title + '.fit'), 1, /SILENT)
    if Keyword_set(Above_ecliptic) then Above_ecliptic_Image_R  = MRDFITS(strcompress(directory + Output_Title + '_Above_ecliptic.fit'), 0, /SILENT)
    if Keyword_set(Above_ecliptic) then Above_ecliptic_Image_CD = MRDFITS(strcompress(directory + Output_Title + '_Above_ecliptic.fit'), 1, /SILENT)
    
    s = size(stackedModel_Image_R)
    x_center=center_in_frame[0]*Output_Size_In_Pixels[0] & y_center=center_in_frame[1]*Output_Size_In_Pixels[1] & z_center=center_in_frame[2]*Output_Size_In_Pixels[2]  
    
  ; Determine the image platescale.  
    ;platescale = float(float(Output_Size_In_Pixels[0]) / (float(plot_range) * body_radius)) ;plate scale of the output image in PIXELS per KM

  ; Compute the projected distance anti-sunward for cases where one object viewed from the center of another object
    if (viewpoint ne 'Above Ecliptic' and viewpoint ne 'Moon Spot') then begin
      cspice_spkpos, body, ephemeris_time, 'J2000', 'NONE', '10', body_position, light_time ; bodys J2000 position with respect to the Sun
      cspice_spkpos, viewpoint, ephemeris_time, 'J2000', 'NONE', '10', observer_position, light_time ; Observers J2000 position with respect to the Sun
      observer_to_body          = body_position-observer_position ; A vector from the observer location to the body
      S_P_O                     = cspice_vsep(body_position, observer_to_body) ;the Sun-body-Observer (STO) angle  
      actual_over_apparent      = 1./(sin(!pi-S_P_O))             ; The projection effect: down tail length is longer than that projected on the sky     
      radii_per_pixel_crosstail = FOV / (S[1]*body_ang_radius) 
      radii_per_pixel_downtail  = FOV / (S[2]*body_ang_radius) 
      pix_per_radii_crosstail   = (1./radii_per_pixel_crosstail) ;(rebinned) to 512x512
      pix_per_radii_downtail    = (1./radii_per_pixel_downtail) / actual_over_apparent
    endif
    
    ; For space-based observers, we'll also need the instrument's boresight direction w/ respect to planet center
      if viewpoint eq 'MESSENGER' then begin
        SC     =   -236     ; MESSENGER spacecraft NAIF ID
        INST   =   -236000  ; MASCS Instrument NAIF ID
        UVVS   =   -236600  ; UVVS atmospheric slit instrument channel MASCS NAIF ID
  
      ; get the FOV and frame definition for the instrument
        cspice_getfov, UVVS, 4, shape, frame, bsight, bounds                                  ; Returns bounds in radians UVVS Atmosphere FOV =  1.0 x 0.04 degrees
        N_bounds = size(bounds, /dim)
  
      ; Convert ephemeris time to spacecraft clock
        cspice_sce2s, SC, ephemeris_time, SC_Clock
  
      ; cspice_ckgp will rbot equire encoded spacecraft clock time, get that encoded time here
        cspice_scencd, SC, SC_Clock, sclkdp
  
      ; Retrieve the 'IAU_Mercury' reference frame to 'INST' reference frame transformation matrix
      ; at time sclkdp with a tolerance of 1.e3 spacecraft clock ticks.
        cspice_ckgp, inst, sclkdp, 1.e3, 'IAU_Mercury', cmat, clkout, found
        if found then begin
          cspice_mtxv, cmat, bsight, boresight_vector_IAU_Mercury                             ; UVVS Boresight Unit vector (IAU Mercury Frame)
          bounds_vectors_IAU_Mercury = dblarr(N_bounds)
          for i = 0, n_bounds[1]-1 do begin
            cspice_mtxv, cmat, bounds[*,i], bounds_vector_IAU_Mercury                         
            bounds_vectors_IAU_Mercury[*,i] = bounds_vector_IAU_Mercury                       ; Unit vectors that bound the aperture's footprint (IAU Mercury Frame)
          endfor
        endif else print, 'No instrument pointing data found' 

      ; Get the planetary coordinates of the tangent point of the boresight ray
        cspice_spkpos, 'MESSENGER', ephemeris_time, 'IAU_Mercury', 'None', body, ptarg, ltime ; ptarg is the Mercury to MESSENGER vector (IAU Mercury Frame)
        SC_to_planet_unit_vector = -ptarg/norm(-ptarg)
        theta = cspice_vsep(boresight_vector_IAU_Mercury, SC_to_planet_unit_vector)           ; angle between the UVVS boresight unit vector and the MESSENGER-planet unit vectors
        SC_to_tangent_point = boresight_vector_IAU_Mercury * norm(-ptarg) * cos(theta)
        tangent_point_MERC_IAU = ptarg + SC_to_tangent_point

      ; And also get the IAU_Mercury planetary coordinates of the corners that bound instrument's aperture
        tangent_point_Bounds_MERC_IAU = dblarr(N_bounds)
        for i = 0, n_bounds[1]-1 do begin
          theta                              = cspice_vsep(bounds_vectors_IAU_Mercury[*,i], SC_to_planet_unit_vector) 
          tangent_point_Bounds_MERC_IAU[*,i] = ptarg + bounds_vectors_IAU_Mercury[*,i] * norm(-ptarg) * cos(theta)
        endfor

        debug = 1 ; Set for inspecting the Planetographic coordinates of the point along the boresight ray that is tangent to the surface
        if Keyword_Set(debug) then begin
          cspice_et2utc, ephemeris_time, 'C', 0, utcstr
          Print, 'At '+utcstr+':' 
          cspice_reclat, tangent_point_MERC_IAU, radius, lon, lat                           ; Validated this result against the MESSENGER PDS DDR
          print, 'UVVS Boresight Tangent Altitude:           ', radius[0] - Body_radius
          print, 'UVVS Boresight Tangent Planetocentric Longitude:', lon[0]*!radeg            
          print, 'UVVS Boresight Tangent Planetocentric Latitude: ', lat[0]*!radeg
          Print, 'Note, however that this simulation uses *PLANETOGRAPHIC* longitude, i.e. with positive to the West.
          Print, 'The PDS DDR give EAST longitude (planetocentric coodinates)'              ; See https://pds-atmospheres.nmsu.edu/PDS/data/messmas_2001/document/uvvs_cdr_ddr_sis.pdf
          cspice_bodvrd, 'Mercury', 'RADII', 3, radii
          flat = (radii[0] - radii[2])/radii[0]
          cspice_recpgr, 'Mercury', tangent_point_MERC_IAU, radii[0], flat, lon, lat, alt   ; Planetographic longitude definition
          print, 'UVVS Boresight Tangent Planetographic Longitude:', lon[0]*!radeg                 
        endif

      ; convert this boresight vector from the IAU Mercury frame in to a body-centered J2000 frame
        cspice_pxform, 'IAU_Mercury', 'J2000', ephemeris_time, IAU_Mercury_to_J2000_transform_matrix
        TP_J2000        = TRANSPOSE(IAU_Mercury_to_J2000_transform_matrix) # tangent_point_MERC_IAU          ; Boresight's tangent point in J2000
        Bounds_TP_J2000 = TRANSPOSE(IAU_Mercury_to_J2000_transform_matrix) # tangent_point_Bounds_MERC_IAU   ; UVVS aperture corners tangent points in J2000

      ; And now rotate the boresight's tangent point from the J2000 frame into the "loc_prime" frame (See Part 1 comments for how this is defined)
        cspice_spkpos, 'Sun', ephemeris_time, 'J2000', 'NONE', body, Body2Sun_vector, light_time
        cspice_spkpos, body, ephemeris_time, 'J2000', 'NONE', viewpoint, Observer2Body_vector, light_time
        cspice_twovec, Observer2Body_vector, 3, Body2Sun_vector, 1, alignment_matrix                         ; Re-defines Observer2Body_vector as 'Z' axis
        TP_prime        = TRANSPOSE(alignment_matrix) # TP_J2000
        Bounds_TP_prime = TRANSPOSE(alignment_matrix) # Bounds_TP_J2000           
        
      ; Move this to a MESSENGER origin, determine its "ra & dec" WRT to the Z axis, then convert this to X & Y image coordinates
        cspice_recrad, TP_prime + [0,0,Body_distance_WRT_obs], TP_prime_dist, TP_prime_RA, TP_prime_dec
        angular_radius = abs(TP_prime_dec - !pi/2.)
        TP_RA  = angular_radius * cos(TP_prime_RA) * 3600.d/CSPICE_RPD()                                     ; This is just an X, Y coordinate in arcseconds from body center
        TP_Dec = angular_radius * sin(TP_prime_RA) * 3600.d/CSPICE_RPD()
        TP_X   = -TP_RA * ang_platescale + x_center                                                          ; RA increases **EASTWARD** not westward in the sky
        TP_Y   = TP_Dec * ang_platescale + y_center
        
      ; And do the same to determine the X and Y pixel coordinates of the boresight bounds   
        Aperture_Corners = dblarr(2, n_bounds[1])
        for i = 0, n_bounds[1]-1 do begin
          cspice_recrad, Bounds_TP_prime[*,i] + [0,0,Body_distance_WRT_obs], Bounds_TP_prime_dist, Bounds_TP_prime_RA, Bounds_TP_prime_dec  
          angular_radius   = abs(Bounds_TP_prime_dec - !pi/2.)
          Bounds_TP_RA     = angular_radius * cos(Bounds_TP_prime_RA) * 3600.d/CSPICE_RPD() 
          Bounds_TP_Dec    = angular_radius * sin(Bounds_TP_prime_RA) * 3600.d/CSPICE_RPD()
          Aperture_Corners[*,i] = [-Bounds_TP_RA * ang_platescale + x_center, Bounds_TP_Dec * ang_platescale + y_center] ; [X, Y] Pixel coordinate of each corner
        endfor

      ; The final boresight pixel. Hack! The physical size and true alignement of the UVVS aperture are neglected here, and including this would give a more fair comparision  
        boresight_pixel = round([TP_X,TP_Y])
    endif  
  
  ; Draw the outline of the body and a line at the terminator for this phase angle
    if keyword_set(label_phase) and viewpoint ne 'Moon Spot' then begin
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

  ; Make an axis in RA and Dec coordinates
    if viewpoint eq 'Moon Spot' then begin
      ;declabs = Jup_dec*!radeg + (FOV/(3600.*N_ticks)) * (findgen(N_ticks + 1) - N_ticks/2.)
      ;RAlabs = Jup_RA*!radeg/15. + (FOV/(3600.*N_ticks)) * (findgen(N_ticks + 1) - N_ticks/2.) * cos(Jup_dec) / 15.
      declabs = anti_Sun_dec*!radeg + (FOV/(3600.*N_ticks)) * (findgen(N_ticks + 1) - N_ticks/2.)
      RAlabs = anti_Sun_RA*!radeg/15. + (FOV/(3600.*N_ticks)) * (findgen(N_ticks + 1) - N_ticks/2.)*cos(anti_Sun_dec) / 15.
      declabs = string(declabs, Format = '(F10.2)')
      RAlabs = reverse(string(RAlabs, Format = '(F10.2)')) ;Sky coordinates: RA Decreases to the left.
    endif else begin
      TICLABELS, Body_dec_WRT_obs*!radeg - FOV/(2.*3600.), N_ticks+1, fov/(60.*(N_ticks)), declabs, delta = 1
      TICLABELS, Body_RA_WRT_obs*!radeg - (FOV*cos(Body_dec_WRT_obs)/(2.*3600.)), N_ticks+1., (fov*cos(Body_dec_WRT_obs)/(60.*N_ticks)) / 15., RAlabs, /RA, delta = 1
    endelse
    RA_axis_format = {XTicks:N_elements(RALABS)-1, XTickname:RALABS, Xtitle:'Right Ascension', $
                      Yticks:N_elements(DecLABS)-1, Ytickname:DecLABS, Ytitle:'Declination', CHARSIZE:1.2, TICKLEN:-0.01}  

  ; Make an axis in Body radii coordinates  
    if viewpoint ne 'Moon Spot' then begin
      xticks = 5 ;xticks = floor(FOV / (tickstep*body_ang_radius))  ;Home many ticks fit in the FOV? 
      xzero = center_in_frame[0]*s[1]
      xtickv = ((tickstep*(findgen(xticks) - (xticks/2)))*pix_per_radii_downtail) + xzero
      xtickname = tickstep*(findgen(xticks) - (xticks/2))
      xtickname = strcompress(string(xtickname, format = '(F10.0)'),/remove_all)  
       
      yticks = 5 ;yticks = floor(FOV / (tickstep*body_ang_radius))
      yzero = center_in_frame[1]*s[2]
      ytickv = ((tickstep*(findgen(yticks) - (yticks/2)))*pix_per_radii_crosstail)+yzero
      ytickname = tickstep*(findgen(yticks) - (yticks/2))
      ytickname = strcompress(string(ytickname, format = '(F10.0)'),/remove_all)  
      
      Body_axis_format = {Xticks:Xticks-1, XTickv:XTickv, XTickname:xtickname, Xtitle:body+' Radii', XSTYLE:1, $
                          Yticks:Yticks-1, YTickv:YTickv, Ytickname:ytickname, Ytitle:body+' Radii', CHARSIZE:1.2, TICKLEN:-0.005}
    endif

  ; Make an axis for the Above Ecliptic viewpoint, body radii coordinates for now    
    if Keyword_set(Above_ecliptic) then begin
      xticks = plot_range/tickstep + 1  
      xzero = center_in_frame[0]*s[1]
      xtickv = tickstep*(findgen(plot_range/tickstep + 1) - center_in_frame[0]*plot_range/tickstep )  *  (Output_size_in_pixels[0]/plot_range) + xzero
      xtickname = tickstep * (findgen(plot_range/tickstep + 1) - center_in_frame[0]*plot_range/tickstep )
      xtickname = strcompress(string(xtickname, format = '(I6)'),/remove_all)               
      
      Yticks = plot_range/tickstep + 1 
      Yzero = center_in_frame[1]*s[2]
      ytickv = tickstep*(findgen(plot_range/tickstep + 1) - .5*plot_range/tickstep )  *  (Output_size_in_pixels[1]/plot_range) + yzero
      ytickname = tickstep * (findgen(plot_range/tickstep + 1) - .5*plot_range/tickstep )
      ytickname = strcompress(string(ytickname, format = '(I6)'),/remove_all)    
      Above_Ecliptic_axis_format = {Xticks:Xticks-1, XTickv:XTickv, XTickname:xtickname, Xtitle:body+' Radii', XSTYLE:1, $
                                    Yticks:Yticks-1, YTickv:YTickv, Ytickname:Ytickname, Ytitle:body+' Radii', CHARSIZE:1.2, TICKLEN:-0.005}                           
    endif                              

  ; Make titles for the plots
    cspice_et2utc, ephemeris_time, "C", 0, utc_current    ; Get a calendar date string
    if viewpoint ne 'Moon Spot' then Title = strcompress(Body + ' - ' + strmid(utc_current,0,17)) $    
    else Title = strcompress( 'Moon Spot - ' + Observatory + ' - '+ strmid(utc_current,0,17))
    Above_Ecliptic_title = strcompress('Projected View From Above Ecliptic - ' + strmid(utc_current,0,17))
     
  ; Label some features in the Moon Spot images
    if viewpoint eq 'Moon Spot' then begin
        Dec_labels = anti_Sun_dec*!radeg + (FOV/(3600.*N_ticks)) * (findgen(N_ticks + 1) - N_ticks/2.)
        RA_labels  = anti_Sun_RA*!radeg/15. + (FOV/(3600.*N_ticks)) * (findgen(N_ticks + 1) - N_ticks/2.)*cos(Anti_sun_dec) / 15.
        ;Dec_labels = Jup_dec*!radeg + (FOV/(3600.*N_ticks)) * (findgen(N_ticks + 1) - N_ticks/2.)
        ;RA_labels  = Jup_RA*!radeg/15. + (FOV/(3600.*N_ticks)) * (findgen(N_ticks + 1) - N_ticks/2.)*cos(Jup_dec) / 15.
        RA_labels  = reverse(RA_labels) ;Eastward increasing
      
      ; Plot Anti-Sun  
        cspice_spkpos, 'Sun', ephemeris_time, 'J2000', 'LT', 'Earth', Sun_Earth_vector, light_time 
        cspice_recrad, -(Sun_Earth_vector-obs_J2000), radius, AntiSun_RA, AntiSun_dec
        AntiSun_pixel = [interpol(findgen(N_ticks + 1) * output_size_in_Pixels[0]/n_ticks,  RA_labels, AntiSun_RA*!radeg/15.), $
                         interpol(findgen(N_ticks + 1) * output_size_in_Pixels[1]/n_ticks, Dec_labels, AntiSun_dec*!radeg)     ]
        print, 'AntiSun_RA, AntiSun_dec', AntiSun_RA*!radeg/15., AntiSun_dec*!radeg                  
    
      ; Plot the Ecliptic Plane 
        cspice_spkpos, 'Sun', ephemeris_time + 30.*24.*(findgen(3600.) - 1800.), 'J2000', 'LT', 'Earth', Sun_Earth_vector, light_time 
        cspice_recrad, -(Sun_Earth_vector-rebin(obs_J2000, 3, 3600.)), radius, Ecliptic_1_RA, Ecliptic_1_dec
        Ecliptic_pixels = [[interpol(findgen(N_ticks + 1) * output_size_in_Pixels[0]/n_ticks,  RA_labels, Ecliptic_1_RA*!radeg/15.)], $
                           [interpol(findgen(N_ticks + 1) * output_size_in_Pixels[1]/n_ticks, Dec_labels, Ecliptic_1_dec*!radeg)  ]]
  
      ; Plot Anti-Moon  
        cspice_spkpos, 'Moon', ephemeris_time, 'J2000', 'LT', 'Earth', Moon_Earth_vector, light_time 
        cspice_recrad, -(Moon_Earth_vector-obs_J2000), radius, AntiMoon_RA, AntiMoon_dec
        AntiMoon_pixel = [interpol(findgen(N_ticks + 1) * output_size_in_Pixels[0]/n_ticks,  RA_labels, AntiMoon_RA*!radeg/15.), $
                          interpol(findgen(N_ticks + 1) * output_size_in_Pixels[1]/n_ticks, Dec_labels, AntiMoon_dec*!radeg)     ]
        print, 'AntiMoon_RA, AntiMoon_dec', AntiMoon_RA*!radeg/15., AntiMoon_dec*!radeg    
      ; Plot Jupiter  
        cspice_spkpos, 'Jupiter', ephemeris_time, 'J2000', 'LT', 'Earth', Jupiter_Earth_vector, light_time 
        cspice_recrad, (Jupiter_Earth_vector-obs_J2000), radius, Jupiter_RA, Jupiter_dec
        Jupiter_pixel = [interpol(findgen(N_ticks + 1) * output_size_in_Pixels[0]/n_ticks,  RA_labels, Jupiter_RA*!radeg/15.), $
                         interpol(findgen(N_ticks + 1) * output_size_in_Pixels[1]/n_ticks, Dec_labels, Jupiter_dec*!radeg)     ]     
        print, 'Jup_RA, Jup_dec', Jup_RA*!radeg/15., Jup_dec*!radeg    
      ; Plot Star 
        Star_RA = ten(10,08,22.31099) & Star_Dec = ten(11,58,01.9516)
        Star_pixel = [interpol(findgen(N_ticks + 1) * output_size_in_Pixels[0]/n_ticks,  RA_labels, Star_RA*!radeg/15.), $
                      interpol(findgen(N_ticks + 1) * output_size_in_Pixels[1]/n_ticks, Dec_labels, Star_dec*!radeg)     ]                      
    endif

;; Write a postscript of the final result, optionally include the body.   
;  cgPS_open, filename = strcompress(directory + Output_title + '_Column_Density.eps'), /ENCAPSULATED, /NOMATCH
;   !P.font=1
;   !p.charsize = 2. 
;   !p.charthick = 2
;   device, SET_FONT = 'Helvetica Bold', /TT_FONT
;   

;   ;Warn if optically thick 
;      optically_thick = where(stackedModel_Image_CD gt Line_data.unity_optical_depth, count)
;      print, count, ' Pixels are optically thick'   
;   ;Scale the column density to be plotted over 1-1000 range
;      Emission_Scaling = strmid(Image_type, 10) 
;      if Emission_Scaling eq 1.e6 then scale = '1.e11'
;      if Emission_Scaling eq 1.e5 then scale = '1.e10'
;      if Emission_Scaling eq 1.e4 then scale = '1.e9'
;      if Emission_Scaling eq 1.e3 then scale = '1.e8'
;      if Emission_Scaling eq 1.e2 then scale = '1.e7'
;      if Emission_Scaling eq 1.e1 then scale = '1.e6'
;      if Emission_Scaling eq 1.e0 then scale = '1.e5'
;      Scale_String = strcompress('x 10!U' + strmid(scale, stregex(scale, 'e') + 1) + '!N')               
;      Scaled_Model_Image_CD  = stackedModel_Image_CD  / float(scale)   
    
;      MinValue = Floor(Min(Scaled_Model_Image_CD))
;      MaxValue = 1.e3 ;Ceil(Max(Scaled_Model_Image_CD))
; 
;      if not keyword_set(Log_Display) then begin
;        loadct, 0, /silent
;        cgImage, Scaled_Model_Image_CD, Stretch=1, MinValue = MinValue, $
;          MaxValue = MaxValue, Position = [0.125, 0.125, 0.9, 0.800], /keep_aspect, $
;          oposition = test, Title = strcompress('Model result - ' + strmid(utc_current,0,17)), charsize =1.2, /Window, /save
;        AXIS,0,0,0, YAXIS = 0, YSTYLE = 1, color=0., charsize=1.3,ythick=3.,$
;            yticks = yticks, ytickv = ytickv, ytickname = ytickname, yticklen = -.01, $
;            Ytitle = STRCOMPRESS(Body + ' Radii'), /data 
;        AXIS,0,0,0, XAXIS = 0, XSTYLE = 1, color=0., charsize=1.3,Xthick=3.,$
;            xticks = xticks, xtickv = xtickv, xtickname = transpose(xtickname), xticklen = -.01, /data, $
;            Xtitle = strcompress('Sunward Distance ' + body + ' Radii')
;        cgColorbar, /Vertical, Position = [!x.window[1]+.07, !y.window[0], !x.window[1]+.1, !y.window[1]], Range=[MinValue, MaxValue], $
;          Title = strcompress('!3' + Particle_data.name + ' Column Density ('+ scale_string +' cm!U-2!N)'), TLocation='Left' ,/right             
;      endif else begin   ; Log scaled image 
;        loadct, 0, /SILENT
;        number_of_colors = 8                                                    ;number of colors / labels on the color bar, loadgood has 8
;        colorbar_label = reverse(MaxValue*.5^findgen(number_of_colors + 1.))    ;color bar labels
;        colorbar_label[0] = 0                                                   ;force the bottom label to zero, true zero has a negaitve inf. scaled value 
;        colorbar_label = strtrim(string(colorbar_label, format = '(i)'), 2)     ;display color bar labels as integers
;        ;Scale the data so that factors of 2 represents integer increases in brightness
;        Log_scale_img_CD = alog(float(Scaled_Model_Image_CD) / MaxValue) / alog(2.) + number_of_colors
;        ;number_of_colors is the value of the new image at the MaxValue peak level
;        cgImage, Log_scale_img_CD, Stretch=1, MinValue = MinValue, MaxValue = number_of_colors, Position = [0.125, 0.125, 0.9, 0.800], /keep_aspect, $
;          oposition = test, Title = strcompress('Model result - ' + strmid(utc_current,0,17)), charsize =1.2, /Window, /save
;        !p.charsize = 2. 
;        AXIS,0,0,0, YAXIS = 0, YSTYLE = 1, color=0., charsize=1.3, ythick=3.,$
;            yticks = yticks, ytickv = ytickv, ytickname = ytickname, yticklen = -.01, $
;            Ytitle = STRCOMPRESS(Body + ' Radii'), /data 
;        AXIS,0,0,0, XAXIS = 0, XSTYLE = 1, color=0., charsize=1.3, Xthick=3.,$
;            xticks = xticks, xtickv = xtickv, xtickname = transpose(xtickname), xticklen = -.01, /data, $
;            Xtitle = strcompress('Sunward Distance ' + body + ' Radii')
;        cgColorbar, /Vertical, Position = [!x.window[1]+.07, !y.window[0], !x.window[1]+.1, !y.window[1]], Range=[MinValue, MaxValue], $
;          Title = strcompress('!3' + Particle_data.name + ' Column Density ('+ scale_string +' cm!U-2!N)'), TLocation='Left', /right, $  
;          ticknames = colorbar_label, divisions = number_of_colors, ticklen = 0, charsize = 1.5
;      
;        if keyword_set(Label_Time) then begin   
;            xyouts, .66,.75, strcompress(Label_Time + ' Hours'), /normal, color=cgColor('white'), charsize = 1.5           
;        endif     
;      endelse           
;  cgPS_Close 

  Log_Display = 1
  cgPS_Open, filename=strcompress(directory + Output_title + '_Emission.eps'), /ENCAPSULATED, /NOMATCH
    !P.font=1
    !p.charsize=1.5
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
      Emission_Scaling = float(strmid(Image_type, 10)) 
      if Emission_Scaling eq 1.e6 then Emission_Label = 'MR'
      if Emission_Scaling eq 1.e5 then Emission_Label = 'x100 kR'
      if Emission_Scaling eq 1.e4 then Emission_Label = 'x10 kR'
      if Emission_Scaling eq 1.e3 then Emission_Label = 'kR'
      if Emission_Scaling eq 1.e2 then Emission_Label = 'x100 R'
      if Emission_Scaling eq 1.e1 then Emission_Label = 'x10 R'
      if Emission_Scaling eq 1.e0 then Emission_Label = 'R'        
      Scaled_Model_Image_R = stackedModel_Image_R / Emission_Scaling ;convert to the Model_Image_R to match the color bar scale  

      minValue = Floor(Min(Scaled_Model_Image_R))
      maxValue = Ceil(Max(Scaled_Model_Image_R))
                  
      if not keyword_set(Log_Display) then begin
        cgImage, Scaled_Model_Image_R, Stretch=1, MinValue = Floor(Min(Scaled_Model_Image_R)), AXKEYWORDS=RA_axis_format, /Axes, $
          MaxValue = Ceil(Max(Scaled_Model_Image_R)), Position = [0.125, 0.125, 0.9, 0.800], /keep_aspect, $
          Title = Title, charsize =1.2, /Window, /save, /noerase
        ;        cgImage, Scaled_Model_Image_R, Stretch=1, MinValue = Min(Scaled_Model_Image_R), AXKEYWORDS=RA_axis_format, /Axes, $
        ;          MaxValue = Max(Scaled_Model_Image_R), Position = [0.125, 0.125, 0.9, 0.800], /keep_aspect, $
        ;          Title = Title, charsize =1.2, /Window, /save, /noerase
        if viewpoint eq 'Moon Spot' then begin
          cgoplot, AntiSun_pixel[0], AntiSun_pixel[1], psym=36, color = 'Yellow'
          cgoplot, AntiMoon_pixel[0], AntiMoon_pixel[1], psym=35, color = 'White'
          cgoplot, Jupiter_pixel[0], Jupiter_pixel[1], psym=16, color = 'Brown'
          cgoplot, Star_pixel[0], Star_pixel[1], psym=46, color = 'White'
          cgoplot, Ecliptic_pixels[*,0], Ecliptic_pixels[*,1], linestyle = 1, color = 'Yellow'  
        endif    
        cgColorbar, /Vertical, Position = [!x.window[1]+.07, !y.window[0], !x.window[1]+.1, !y.window[1]], Range=[MinValue, MaxValue], $
          Title = Strcompress('Simulated ' + Line_data.name + ' Line Emission ('+Emission_Label+')'), TLocation='Left', /right, charsize = 1.5, /noerase, TICKINTERVAL = maxValue/8.              
      endif else begin   ; Log scaled image
        number_of_colors = 8                                                    ;number of colors / labels on the color bar, loadgood has 8
        colorbar_label = reverse(MaxValue*.5^findgen(number_of_colors + 1.))    ;color bar labels
        colorbar_label[0] = 0                                                   ;force the bottom label to zero, true zero has a negaitve inf. scaled value 
        colorbar_label = strtrim(string(colorbar_label, format = '(i)'), 2)     ;display color bar labels as integers
        ;Scale the data so that factors of 2 represents integer increases in brightness
        Log_scale_img_R = alog(float(Scaled_Model_Image_R) / MaxValue) / alog(2.) + number_of_colors
        ;number_of_colors is the value of the new image at the MaxValue peak level
        cgImage, Log_scale_img_R, Stretch=1, MinValue=MinValue, MaxValue = number_of_colors, Position = [0.075, 0.125, 0.9, 0.90], AXKEYWORDS=RA_axis_format, /Axes, /KEEP_ASPECT, $
          Title = Title, charsize =1.2, /Window, /save
        cgColorbar, /Vertical, Position = [!x.window[1]+.07, !y.window[0], !x.window[1]+.1, !y.window[1]], Range=[MinValue, MaxValue], $
          Title = Strcompress('Simulated ' + Line_data.name + ' Line Emission ('+Emission_Label+')'), TLocation='Left', /right, $
          ticknames = colorbar_label, divisions = number_of_colors, ticklen = 0, charsize = 1.5 
      endelse
  cgPS_Close    

  ;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  if Keyword_set(Above_ecliptic) then begin
    cgPS_Open, filename=strcompress(directory + Output_title + '_Above_ecliptic.eps'), /ENCAPSULATED, /NOMATCH
    !P.font=1
    !p.charsize=1.5
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
      Emission_Scaling = strmid(Image_type, 10) 
      if Emission_Scaling eq 1.e6 then Emission_Label = 'MR'
      if Emission_Scaling eq 1.e5 then Emission_Label = 'x100 kR'
      if Emission_Scaling eq 1.e4 then Emission_Label = 'x10 kR'
      if Emission_Scaling eq 1.e3 then Emission_Label = 'kR'
      if Emission_Scaling eq 1.e2 then Emission_Label = 'x100 R'
      if Emission_Scaling eq 1.e1 then Emission_Label = 'x10 R'
      if Emission_Scaling eq 1.e0 then Emission_Label = 'R'        
      Above_ecliptic_Image_R = Above_ecliptic_Image_R / Emission_Scaling ;convert to the Model_Image_R  to match the color bar scale  
        
      minValue = Floor(Min(Above_ecliptic_Image_R))
      maxValue = Ceil(Max(Above_ecliptic_Image_R))
                  
      if not keyword_set(Log_Display) then begin
        cgImage, Above_ecliptic_Image_R, Stretch=1, MinValue = Floor(Min(Above_ecliptic_Image_R)), $
          MaxValue = Ceil(Max(Above_ecliptic_Image_R)), Position = [0.125, 0.125, 0.9, 0.800], AXKEYWORDS=Above_Ecliptic_axis_format, /Axes, /KEEP_ASPECT_RATIO, $
          oposition = test, Title = Above_Ecliptic_title, charsize =1.2, /Window, /save
        cgColorbar, /Vertical, Position = [!x.window[1]+.07, !y.window[0], !x.window[1]+.1, !y.window[1]], Range=[MinValue, MaxValue], $
          Title = Strcompress('Simulated ' + Line_data.name + ' Line Emission ('+Emission_Label+')'), TLocation='Left', /right, charsize = 1.5              
      endif else begin   ; Log scaled image
        number_of_colors = 8                                                    ;number of colors / labels on the color bar, loadgood has 8
        colorbar_label = reverse(MaxValue*.5^findgen(number_of_colors + 1.))    ;color bar labels
        colorbar_label[0] = 0                                              ;force the bottom label to zero, true zero has a negaitve inf. scaled value 
        ;colorbar_label = strtrim(string(colorbar_label, format = '(i)'), 2)     ;display color bar labels as integers
        colorbar_label = strcompress(string(colorbar_label, format = '(F5.2)'), /remove_all)
        ;Scale the data so that factors of 2 represents integer increases in brightness
        Log_scale_img_R = alog(float(Above_ecliptic_Image_R) / MaxValue) / alog(2.) + number_of_colors
        ;number_of_colors is the value of the new image at the MaxValue peak level
        cgImage, Log_scale_img_R, Stretch=1, MinValue=MinValue, MaxValue = number_of_colors, Position = [0.125, 0.125, 0.9, 0.800], AXKEYWORDS=Above_Ecliptic_axis_format, /Axes, /KEEP_ASPECT_RATIO, $
          Title = Above_Ecliptic_title, charsize =1.2, /Window, /save
        cgColorbar, /Vertical, Position = [!x.window[1]+.07, !y.window[0], !x.window[1]+.1, !y.window[1]], Range=[MinValue, MaxValue], $
          Title = Strcompress('Simulated ' + Line_data.name + ' Line Emission ('+Emission_Label+')'), TLocation='Left', /right, $
          ticknames = colorbar_label, divisions = number_of_colors, ticklen = 0, charsize = 1.5 
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