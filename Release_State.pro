pro Release_State, loc, speed_distribution, surface_distribution, speed

COMMON Model_shared, Body, Ephemeris_time, Seed, Directory, Particle_data, Line_data, Debug
COMMON Output_shared, Plot_range, Output_Size_In_Pixels, Output_Title, Center_in_frame, viewpoint, FOV, N_ticks, Tickstep, Observatory, Above_Ecliptic, Boresight_Pixel, Aperture_Corners

  N_particles = N_elements(loc[0,*])

  ; Calculate the flatness coefficient and planetary radius
    cspice_bodvrd, body, 'RADII', 3, radii
    re   = radii[0]
    rp   = radii[2]
    flat = ( re - rp ) / re

; ================================================= LAUNCH SPEEDS =================================================================================================

; Parse up the Main level input string get the distribution of speeds the particles will launch with:
  if STRPOS(speed_distribution, 'Step') eq 0 then begin
    V_start = STRPOS(speed_distribution, '[') 
    comma = STRPOS(speed_distribution, ',') 
    V_end = STRPOS(speed_distribution, ']') 
    Vmin = float(STRMID(speed_distribution, V_start+1, comma-v_start-1))
    Vmax = float(STRMID(speed_distribution, comma+1, v_end-comma-1))
    ;a random distribution of speeds in the specified input range for each particle:
    speed = randomu(seed,N_particles) * (Vmax - Vmin) + Vmin 
    isotropic_release = 1 
  endif
  if STRPOS(speed_distribution, 'Maxwellian') eq 0 then begin
    end_of_temperature_string_position = STRPOS(speed_distribution, 'K') 
    Temperature = float(STRMID(speed_distribution,11,end_of_temperature_string_position-1))
    Distribution = 'Maxwellian'
    distribution_parameters = {Maxwellian_Structure, Distribution:Distribution, Temperature:Temperature}
    speed = generate_velocity_distribution(N_particles, distribution_parameters)
    isotropic_release = 1 ;ejection angles randoimly distributed over 2 pi steradians
  endif
  if STRPOS(speed_distribution, 'MBF') eq 0 then begin
    end_of_temperature_string_position = STRPOS(speed_distribution, 'K') 
    Temperature = float(STRMID(speed_distribution,4,end_of_temperature_string_position-1))
    Distribution = 'MBF'
    distribution_parameters = {Maxwellian_Structure, Distribution:Distribution, Temperature:Temperature}
    speed = generate_velocity_distribution(N_particles, distribution_parameters)
  endif
  if STRPOS(speed_distribution, 'Cone') eq 0 then begin
    end_of_temperature_string_position = STRPOS(speed_distribution, 'K')
    Temperature = float(STRMID(speed_distribution,4,end_of_temperature_string_position-1))
    Temperature = 3500. ;Hack!!! hard code this!
    Distribution = 'MBF'
    distribution_parameters = {Maxwellian_Structure, Distribution:Distribution, Temperature:Temperature}
    speed = generate_velocity_distribution(N_particles, distribution_parameters)
  endif
  if STRPOS(speed_distribution, 'Drifting_MBF') eq 0 then begin
    end_of_temperature_string_position = STRPOS(speed_distribution, 'K') 
    Temperature = float(STRMID(speed_distribution,13,end_of_temperature_string_position-1))
    Distribution = 'MBF'
    distribution_parameters = {Maxwellian_Structure, Distribution:Distribution, Temperature:Temperature}
    speed = generate_velocity_distribution(N_particles, distribution_parameters)
  endif
  if STRPOS(speed_distribution, 'Kappa') eq 0 then begin
    end_of_temperature_string_position = STRPOS(speed_distribution, 'K_') 
    begin_kappa_string_position = STRPOS(speed_distribution, 'K=') 
    Temperature = float(STRMID(speed_distribution,6,end_of_temperature_string_position-1))
    Kappa       = float(STRMID(speed_distribution,begin_kappa_string_position+2))
    Distribution = 'Kappa'
    distribution_parameters = {Kappa_Structure, Distribution:Distribution, Temperature:Temperature, Kappa:Kappa}
    speed = generate_velocity_distribution(N_particles, distribution_parameters)
  endif

; ============================================== LAUNCH COORDINATES =================================================================================================
; Get the particle launch locations in body-centered J2000 coordinates
  CASE 1 OF
    STRMATCH(surface_distribution, 'Global'): BEGIN                                                           ; Generate isotropic latitudes and longitudes for all particles
        lat = asin(1.d - 2.d*randomu(seed, N_particles, /double))             ; gives a random number, -90 to 90 degrees (in radians) weighted towards the equator (more surface area per unit latitude)
        lon = 360.d*randomu(seed, N_particles, /double) /!RADEG               ; gives random longitudes 0 to 360 (in radians)
        cspice_pgrrec, body, lon, lat, replicate(0.D, N_particles), re, flat, release_points_body_fixed; convert these to cartesian, body fixed coords untis of planetary radii here, not km
        release_points_J2000 = release_points_body_fixed / re                 ; since these are random and isotropic in x, y, z, the reference frame does not matter, use body radii coordinates
    END  
    STRMATCH(surface_distribution, 'Dayside'): BEGIN
      release_points_J2000      = dblarr(3, N_Particles)
      release_points_Body_fixed = dblarr(3, N_Particles)
      cspice_subpnt, 'Near point: ellipsoid', body, ephemeris_time, 'IAU_'+body, 'LT', viewpoint, sub_observer_point_planet_frame, trgepc, srfvec      
      cspice_recpgr, body, sub_observer_point_planet_frame, re, flat, sub_observer_lon, sub_observer_lat, sub_observer_radius 
      ss_lon = fltarr(N_Particles) & ss_lat = fltarr(N_Particles) 
      for i = 0, N_particles - 1 do begin ; sub-solar point moves depending on the release time, so we need a for loop to do this... 
        
        ; Get an x,y,z subsolar point (in km from planet center) in the body-fixed frame at the time of release.
          cspice_subslr, 'Near point: ellipsoid', body, ephemeris_time - loc[9,i], 'IAU_'+body, 'none', viewpoint, sub_solar_point_planet_frame, trgepc, srfvec

        ; Convert this sub-solar point to planetographic latitude and longitude, in the 'IAU_' body fixed frame
          cspice_recpgr, body, sub_solar_point_planet_frame, re, flat, sub_solar_lon, sub_solar_lat, sub_solar_point_radius        
          
        ; Generate a randomn latitude and longitude so that the zenith of the hemisphere is along the positive y direction.
          lat = asin(1.0-2.0*randomu(seed))                            ; gives a random number between -90 and 90 degrees (in radians) 
                                                                       ; weighted towards the equator (more surface area per unit latitude)
          lon = ((180.*randomu(seed)))/!radeg - !pi/2. + sub_solar_lon ; gives a random longitude -90 and 90 degrees from the sub-solar point(in radians)

          ss_lat[i] = lat 
          ss_lon[i] = Lon - sub_solar_lon

        ; Convert the release point from planetographic coordinates to rectangular coordinates, still in the body fixed frame
          cspice_pgrrec, body, lon, lat, 0., re, flat, release_point_body_fixed

        ; get the rotation matrix to J2000 frame the rest of the model runs in.
          cspice_pxform, 'IAU_'+body, 'J2000', ephemeris_time - loc[9,i], body_fixed_to_J2000_rotation_matrix
          cspice_mxv, body_fixed_to_J2000_rotation_matrix, release_point_Body_fixed, release_point_J2000_frame
          release_points_J2000[*,i]      = release_point_J2000_frame / re  ; use body radii as the coodinates here
          release_points_body_fixed[*,i] = release_point_Body_fixed
      endfor

      ; normalize for the photon flux per unit surface area
        SZA_map = cos(SS_Lon)*cos(SS_Lat)
        SZA_map[where(SZA_map le 0., /NULL)] = 0.
        normalize_for_SZA = float(n_particles) / total(SZA_map)  
        loc[4,*] = normalize_for_SZA * SZA_map
        
    END
    STRMATCH(surface_distribution, 'Point*'): BEGIN
       
        release_points_J2000      = dblarr(3, N_Particles)
        release_points_Body_fixed = dblarr(3, N_Particles) 
         
;      ; Get the sub-observer coordinates
;        cspice_subpnt, 'Near point: ellipsoid', body, ephemeris_time, 'IAU_'+body, 'LT', viewpoint, sub_observer_point_planet_frame, trgepc, srfvec
;        cspice_recpgr, body, sub_observer_point_planet_frame, re, flat, sub_observer_lon, sub_observer_lat, sub_observer_radius ; Convert to planetographic latitude and longitude, in the 'IAU_' body fixed frame
;
;      ; Get an x,y,z subsolar point (in km from planet center) in the body-fixed frame at the time of release.
;        cspice_subslr, 'Near point: ellipsoid', body, ephemeris_time, 'IAU_'+body, 'none', viewpoint, sub_solar_point_planet_frame, trgepc, srfvec ; not sure whay we'd need light time 
;        cspice_recpgr, body, sub_solar_point_planet_frame, re, flat, sub_solar_lon, sub_solar_lat, sub_solar_radius ; Convert to planetographic latitude and longitude, in the 'IAU_' body fixed frame
;        lon  = (sub_solar_lon + !pi/2. + 2.*!pi) mod (2.*!pi)
;        lat  = 0.d * cspice_rpd()

      ; Set the particle release points to a specified location, try equatorial dawn for a meteor impact
      ; Set a planetographic longitude, latitude, altitude position.
        alt  = 0.d
        lon = float(STRMID(surface_distribution, STRPOS(surface_distribution, '[')+1, STRPOS(surface_distribution, ',')-STRPOS(surface_distribution, '[')-1)) * cspice_rpd()
        lat = float(STRMID(surface_distribution, STRPOS(surface_distribution, ',')+1, STRPOS(surface_distribution, ']')-STRPOS(surface_distribution, ',')-1)) * cspice_rpd()

      ; Convert lat, lon, altitude to planetographic rectangular coordinates
        cspice_pgrrec, body, lon, lat, alt, re, flat, release_point_body_fixed
      
      ; Convert planetographic coodinates to the J2000 frame at the ephemeris time for each particles release
        for i = 0, N_Particles - 1 do begin
          cspice_pxform, 'IAU_'+body, 'J2000', ephemeris_time - loc[9,i], frame_rotate_matrix
          cspice_mxv, frame_rotate_matrix, release_point_body_fixed, release_point_J2000_frame
          release_points_J2000[*,i]      = release_point_J2000_frame / radii[0] ;the initial packet starting points in cartesian units of planetary radii
          release_points_body_fixed[*,i] = release_point_body_fixed
        endfor

    END
    STRMATCH(surface_distribution, 'From_Map*'): BEGIN
      map_filename = 'C:\IDL\Generic Model V2\read_write\Surface_reservoir\Mean_Reimpacting_Flux_Map.fit'
      map = mrdfits(map_filename) ; TAKE CARE THAT LONGITUDE IS IN THE PROPER CONVENTION HERE: planetographic WEST longitude is the x axis (left-handed longitude)      
      
      ; Get an x,y,z subsolar point (in km from planet center) in the body-fixed frame at the time of release.
      sub_solar_point_planet_frame = fltarr(3, N_Particles)
      for i = 0, N_Particles - 1 do begin
        cspice_subslr, 'Near point: ellipsoid', body, ephemeris_time - loc[9,i], 'IAU_'+body, 'none', viewpoint, sub_solar_point_BF, trgepc, srfvec
        sub_solar_point_planet_frame[*,i] = sub_solar_point_BF
      endfor
      
      ; Convert this sub-solar point to planetographic latitude and longitude, in the 'IAU_' body fixed frame
        cspice_recpgr, body, sub_solar_point_planet_frame, re, flat, sub_solar_lon, sub_solar_lat, sub_solar_point_radius

      ; Generate a randomn latitude and longitude so that the zenith of the hemisphere is along the positive y direction.
        lat = asin(1.0-2.0*randomu(seed, N_particles, /double))                            ; gives a random number between -90 and 90 degrees (in radians)
        ; weighted towards the equator (more surface area per unit latitude)
        lon = ((180.*randomu(seed, N_particles, /double)))/!radeg - !pi/2. + sub_solar_lon ; gives a random longitude -90 and 90 degrees from the sub-solar point(in radians)
      
        lon = lon mod (2.*!pi)
;      ; Generate isotropic release points in the body-fixed frame
;        lat = asin(1.d - 2.d*randomu(seed, N_particles, /double))             ; gives a random number, -90 to 90 degrees (in radians) weighted towards the equator (more surface area per unit latitude)
;        lon = 360.d*randomu(seed, N_particles, /double) /!RADEG               ; gives random longitudes 0 to 360 (in radians)

      ss_lat = lat
      ss_lon = Lon - sub_solar_lon
      
        cspice_pgrrec, body, lon, lat, replicate(0.D, N_particles), re, flat, release_points_body_fixed; convert these to cartesian, body fixed coords units of planetary radii here, not km

      ; Convert planetographic coodinates to the J2000 frame at the ephemeris time for each particles release
        release_points_J2000      = dblarr(3, N_Particles)
        for i = 0, N_Particles - 1 do begin
          cspice_pxform, 'IAU_'+body, 'J2000', ephemeris_time - loc[9,i], frame_rotate_matrix
          cspice_mxv, frame_rotate_matrix, release_points_body_fixed[*,i], release_point_J2000_frame
          release_points_J2000[*,i]      = release_point_J2000_frame / radii[0] ;the initial packet starting points in cartesian units of planetary radii
        endfor
        
      ;------------------------------------hack------------------------------------------------------------------------------  
        ;normalize = rebin(transpose(mean(map, dim = 1)), 360, 180)           ; HACK this normalization will remove the latitude Na
        normalize = mean(map)
        
        map = map / normalize                                                 ; normalize the map to unity
        map[*,179] = 1.                                                       ; fix the edge
        map[where(map eq 0, /null)] = 1.                                      ; fix any missing pixels
        map[83:97,*]   = 8. * map[83:97,*]                                    ; hacking here
        map[263:277,*] = 8. * map[263:277,*]
      ;------------------------------------hack------------------------------------------------------------------------------    
      loc[4,*] = interpolate(map, lon*!radeg, lat*!radeg + 90.)               ; weight the particle initial contents by the map
 
      ; normalize for the photon flux per unit surface area
        SZA_map = cos(SS_Lon)*cos(SS_Lat)
        SZA_map[where(SZA_map le 0., /NULL)] = 0.
        normalize_for_SZA = float(n_particles) / total(SZA_map)  
        loc[4,*] = loc[4,*]*normalize_for_SZA * SZA_map

    END  
    ELSE: PRINT, 'Undefined spatial distribution requested as an input'  
  ENDCASE

  ; The final coordinates. Units of Body-radii. Body-centered J2000 coordinates
    loc[0:2,*] = temporary(release_points_J2000) ; Allocate the release coordinate to the initialize the loc array

  if keyword_set(debug) then begin
    angles_from_normal = fltarr(N_particles) ; for inspecting the ejection angle wrt local surface normal vector
    cspice_recpgr, body, release_points_body_fixed, radii[0], flat, longitudes, latitudes, radius ; Convert to planetographic latitude and longitude, in the 'IAU_' body fixed frame
    cgplot, (longitudes*!radeg + 360.) mod 360., latitudes*!radeg, psym=3, xr=[0,360], yr=[-90,90], title='Red: Sub-Sol Blue: Sub-Obs' ; particle launch coordinates 
    cgplot, (sub_solar_lon*!radeg + 360.) mod 360., sub_solar_lat*!radeg, psym = 16, symsize = 2, color = 'red', /overplot             ; the sub-solar point
    cgplot, (sub_observer_lon*!radeg + 360.) mod 360., sub_observer_lat*!radeg, psym = 16, symsize = 2, color = 'blue', /overplot      ; the sub-observer poin
  endif

; ============================================== LAUNCH DIRECTIONS =================================================================================================

if keyword_set(isotropic_release) then begin ;Generate randomly distributed isotropic release positions and velocities on a unit sphere 

;Explicity specified launch angles
  ;  Theta = 90. / !radeg ;Theta = zenith angle WRT surface normal, units are radians
  ;  Phi = (90.) / !radeg ;Phi = azimuth angle WRT surface normal, units are radians
  ;  V = fltarr(3, N_particles)
  ;  if keyword_set(debug) then angle_from_surface_normal = fltarr(N_particles)
  ;  for i = 0, N_particles-1 do begin  
  ;    ;step 1: get any vector perpendicular to the surface vector by just taking the cross product with the z axis
  ;      cspice_vcrss, loc[0:2,i], [0,0,1.1], X_Cross_Z 
  ;    ;step 2: spin the cross-product vector about the surface vector by the angle phi
  ;      cspice_vrotv, X_Cross_Z, loc[0:2,i], phi, X_Cross_V_phi_spun
  ;    ;step 3: now that azimuth (phi) has been determined, rotate the surface vector by theta to get the velocity vector
  ;      cspice_vrotv, loc[0:2,i], X_Cross_V_phi_spun, theta, velocity_vector
  ;      V[*,i] = velocity_vector
  ;    ;step 4: if desired, verify the result  
  ;      if keyword_set(debug) then angle_from_surface_normal[i] = !radeg * cspice_vsep(loc[0:2,i], V[*,i])
  ;  endfor
  ;  if keyword_set(debug) then cghistoplot, angle_from_surface_normal, binsize = 1., title = 'Ejection Angle Relative to Surface Normal'
  ;  loc[5:7,*] = temporary(V)

;Random isotropic launch angles upward into 2 Pi sterad from the surface normal    
 
  ; generate random points of unity length for the velocity vectors:
    theta = asin(1.d - 2.d*randomu(seed,N_particles, /double)) ; gives random numbers evenly distributed between -90 and 90 degrees (in radians) 
    phi   = 360.*randomu(seed,N_particles,/double)/!RADEG      ; gives random numbers evenly distributed between 0 and 360 degrees (in radians)
    loc[5,*]=cos(theta)*cos(phi)  ;the x-direction, a random number from -1 to 1
    loc[6,*]=cos(theta)*sin(phi)  ;the y-direction, a random number from -1 to 1
    loc[7,*]=sin(theta)           ;the z-direction, a random number from -1 to 1 
    
    ;Now there are two sets of vectors (position and velocity), both randomly distributed over 4 pi steradians on a unit sphere.
    ;To get the velocity vectors to launch upwards, correlate these two arrays. . .
    ;For an isotropic velocity over distribution 2 pi steradians, flip the velocity vectors which are more than 90 degrees from the surface normal directions
    ;This makes the velocity vectors randomly point skyward instead of planet-ward
    for q = 0, N_particles-1 do begin
      angle = cspice_vsep(loc[0:2,q],loc[5:7,q]) ;for each particle, compute the angle between surface normal and velocity vectors.
      if (angle gt !pi/2.) then loc[5:7,q] = -1. * loc[5:7,q] ;if at all planet-ward, flip it (invert 180 degrees) so that it points skyward.
      if keyword_set(debug) then angles_from_normal[q] = !Radeg * cspice_vsep(loc[0:2,q], loc[5:7,q])
    endfor 
    if keyword_set(debug) then begin ;check to be sure these have a cosine x sine like dependence
      window, 0, title='Ejection Angles Relative to Surface Normal'
      cgHISTOPLOT, angles_from_normal, xtitle = 'Ejection Angle from Surface Normal (Degrees)', Ytitle = 'Number of particles', $
        Histdata = histdata 
      ;To help guide the eye:   
      oplot, !Radeg*(!pi/2.) * (findgen(10000.)/10000.), max(histdata) * sin((!pi/2.) * (findgen(10000.)/10000.)), color = 255.
      ;The ejection angles from the exobase surface normal should have a median value of 60 for an isotropic maxwell boltzmann distribution.  
    endif
endif 

;---------------------------------MAXWELL BOLZMANN FLUX DISTRIBUTIONS ARE *anisotropic* PREFERENTIALLY RADIAL (I.E. NOT ISOTROPIC)--------------------------------------------------------------------------------- 
; weight ejectrion angles as cos of the surface normal
  if STRPOS(speed_distribution, 'MBF') eq 0 then begin
  
    angles_from_normal = fltarr(N_particles)
  
    ; Ejection angle's direction weighted by the cosine of the angle between ejection and the surface normal.
    
      r     = findgen(N_particles) / N_particles                      ; a zero to 1 array of n_particle increments
      theta = fltarr(N_particles)                                     ; the independent variable that we're weighting over, in this case the angle each particle will have 
  
      ; Define the probability distribution function that we're going to weight over.
        distribution_function = sin(!pi*r/2.) * cos(!pi*r/2.)  

      sum = TOTAL( distribution_function, /CUMULATIVE )               ; get the cumulative sum of the distribution over elements 0 to i, ie, the cumulative distribution function
      sum = sum/max(sum)                                              ; normalize the CDF to a maximum of one

      for i = 1, N_particles-1 do begin                               ; step through the CDF, I don't see a way to vectorize this
        inrange = where((r ge sum[i-1]) and (r lt sum[i]), count, /null)    
        if count gt 0 then theta[inrange] = i
      endfor
      theta = (theta/N_particles) * (!pi/2.)                          ; scale the distribution from 0 to pi/2 radians.

      ; -> randomize them so that they differ for each particle following the random number seed.
        randomNumbers = RANDOMU(seed, N_particles) 
        scrambled_indicies = SORT(randomNumbers)
        theta[indgen(N_particles)] = theta[scrambled_indicies]

      phi = 360.*randomu(seed, N_particles) / !RADEG                  ; Isotropic in azimuth, random numbers evenly distributed between 0 and 360 degrees (in radians)
  
      ; We're in a spherical coordinate system where z is defined as polar_angle_distribution (theta) = 0 is the z axis, phi is the angle relative to the x axis.
        loc[5,*] = sin(theta)*cos(phi)  ;the x-direction 
        loc[6,*] = sin(theta)*sin(phi)  ;the y-direction
        loc[7,*] = cos(theta)           ;the z-direction  
        
      ; Okay, now a cosine dependent distribution of polar angles is created as cartesian vectors. 
      ; However this distribution is relative an aribitrary coordinate system, where z (loc[7,*]) is up, not the individual surface normal vectors xyz ->
      ; Compute each surface normals rotation from absolute x,y,z, then apply that rotation to the velocity vectors one by one, 
      ; Making each velocity vector specific to a particular surface location.
      ; IT'D BE USEFUL TO VECTORIZE THIS!! (instead of the for loop below) 
        z_axis = [0., 0., 1.] ; define the z-axis in vector space.
   
        ;Rotate to align the z axis in the velocity vector distribution to the surface normal vectors. . .     
        for q = 0, N_particles-1 do begin
          Surface_Normal              = loc[0:2,q] ; Vectors normal to the local surface.      
          Velocity_Vector_WRT_Surface = loc[5:7,q] ; Velocity unit vectors in an arbitary Z is up frame, we want to rotation these into a Surface normal is up frame
          cspice_vcrss, surface_normal, z_axis, surface_normal_cross_z           ; the cross product of the two. This resultant vector is the axis of rotation 
          Suface_Normal_Angle_WRT_Z = cspice_vsep(surface_normal,z_axis)         ; the angle between z axis and the surface vector where the particle is released.       
          
          ; rotate the velocity vector so that the z axis (which it was generated WRT) is now the surface normal. 
            cspice_vrotv, Velocity_vector_WRT_Surface, surface_normal_cross_z, -(Suface_Normal_Angle_WRT_Z), Velocity_vector_WRT_Absolute 
            loc[5:7,q] = Velocity_vector_WRT_Absolute                            ; put it in the loc array velocity corresponding to particle q           
          if keyword_set(debug) then angles_from_normal[q] = !Radeg * cspice_vsep(loc[0:2,q], loc[5:7,q]) ; for each particle compute the ejection angle relative to normal as a check   
        endfor 

        ;check to be sure these have a cosine like dependence . . . should peak at 45 degrees       
        if keyword_set(debug) then cgHISTOPLOT, angles_from_normal, binsize = 2.,xtitle = 'Ejection Angle from Surface Normal (Degrees)', Ytitle = 'Number of particles' 
  endif
  
;---------------------------------CONE-SHAPED EJECTA ARE KNOWN IN METEOR IMPACT DUST--------------------------------------------------------------------------------- 
; cf. Horanyi et al. 2015 Nature 
;     Bernardoni, E. A., Szalay, J. R., & Horányi, M. (2018). Impact Ejecta Plumes at the Moon. Geophysical Research Letters. doi:10.1029/2018gl079994 
;     Szalay et al. (2018)
if STRPOS(speed_distribution, 'Cone') eq 0 then begin
  
  ;debug = 1
  
  angles_from_normal = fltarr(N_particles)

  ; Ejection direction is weighted.
  ; Get the distribution function we're going to weight over.
    distribution_function = sin((!pi/2.) * findgen(N_particles)/N_particles) * cos((!pi/2.)*findgen(N_particles)/N_particles)^3. ; Hack! 
  ; But looks okay... compare with Szalay et al 2018, Impact Ejecta and Gardening in the Lunar Polar Regions
  
  sum = distribution_function ;sum up the distribution function element by element into an array
  for i = 1, n_elements(distribution_function)-2 do begin
    sum(i) = sum(i-1)+((distribution_function(i-1) + distribution_function(i))/2.)
  endfor
  sum = sum/max(sum) ;normalize to a maximum of one
  distribution_function = distribution_function/float(max(sum))
  r = findgen(N_particles)/N_particles
  theta = fltarr(N_particles) ;speed is the speed each particle will have
  i = long(0)
  while i le n_elements(sum)-2 do begin ;step through the sum of the distribution
    top = sum(i)
    bottom = 0.
    if i gt 0 then bottom = sum(i-1) ;the lower bound of the sum of the distribution
    inrange = where((r lt top) and (r ge bottom))
    if inrange(0) ne -1 then theta(inrange) = i
    i=i+long(1)
    if sum(i-1) eq 1. then i = n_elements(sum)-1 ; when stepping through the distribution totals to 1, stop the loop
  endwhile
  theta = theta*((!pi/2.)/N_particles) ;scale it from 0 to pi/2.

  ; -> randomize them so that they differ for each particle following the random number seed.
  randomNumbers = RANDOMU(seed, N_particles)
  scrambled_indicies = SORT(randomNumbers)
  theta[indgen(N_particles)] = theta[scrambled_indicies]

  phi = 360.*randomu(seed,N_particles)/!RADEG  ;gives random numbers evenly distributed between 0 and 360 degrees (in radians), random and isotropic in azimuth

  ;where in a spherical coordinate system where z is defined as polar_angle_distribution (theta) = 0 is the z axis, phi is the angle relative to the x axis.
  loc[5,*]=sin(theta)*cos(phi)  ;the x-direction
  loc[6,*]=sin(theta)*sin(phi)  ;the y-direction
  loc[7,*]=cos(theta)           ;the z-direction

  ; Okay, now a cosine dependent distribution of polar angles is created as cartesian vectors.
  ; However this distribution of vectors is relative to the individual surface vectors xyz and not the absolute coordinates. ->
  ; Compute each surface normals rotation from absolute x,y,z, then apply that rotation to the velocity vectors one by one,
  ; Making each velocity vector specific to a particular surface location.
  ; IT'D BE USEFUL TO VECTORIZE THIS!! (instead of the for loop below)
  z_axis = [0., 0., 1.] ;define the z-axis in vector space.

  ;Rotate to align the z axis in the velocity vector distribution to the surface normal vectors. . .
  for q = 0, N_particles-1 do begin
    Surface_Normal              = loc[0:2,q] ; Vectors normal to the local surface.
    Velocity_Vector_WRT_Surface = loc[5:7,q] ; ?????????in the surface normal frame before rotation to the absolute.????????
    cspice_vcrss,surface_normal,z_axis,surface_normal_cross_z ;the cross product of the two. This resultant vector is the axis of rotation
    Suface_Normal_Angle_WRT_Z = cspice_vsep(surface_normal,z_axis) ;the angle between z axis and the surface vector where the particle is released.
    ;rotate the velocity vector so that the z axis (which it was generated WRT) is the surface normal.
    cspice_vrotv, Velocity_vector_WRT_Surface, surface_normal_cross_z, -(Suface_Normal_Angle_WRT_Z), Velocity_vector_WRT_Absolute
    loc[5:7,q] = Velocity_vector_WRT_Absolute ;put it in the loc array velocity corresponding to particle q
    if keyword_set(debug) then angles_from_normal[q] = !Radeg * cspice_vsep(loc[0:2,q], loc[5:7,q]) ; for each particle compute the ejection angle relative to normal as a check
  endfor
  ;check to be sure these have a cosine like dependence . . . should peak at 45 degrees
  if keyword_set(debug) then cgHISTOPLOT, angles_from_normal, binsize = 2.,xtitle = 'Ejection Angle from Surface Normal (Degrees)', Ytitle = 'Number of particles'

endif  
return
end