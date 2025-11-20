function RK_4_acceleration, Body_xyz, Sun_xyz, time, v

  ; Body_xyz is a 3 x N_particles vector of the J2000 particle locations with respect to Body center at 'time'
  ; Sun_xyz  is a 3 x N_particles vector of the J2000 particle locations with respect to Sun center at 'time'
  ; time is the particle's time before (emphemeris_time - obs_body_ltime)
  ; v is the heliocentric velocity in km/s
  ; 

  COMMON Gravity, GM_Body, GM_Sun, GM_Parent
  COMMON Model_shared, Body, Ephemeris_time, Obs_Body_Ltime, Parent_ID, Seed, Directory, Particle_data, Line_data, Debug

    x = Body_xyz[0,*]
    y = Body_xyz[1,*]
    z = Body_xyz[2,*]
    sun_particle = Sun_xyz

  ; Get the distance cubed from the Body and the Sun
    r_Body3 = ([x^2] + [y^2] + [z^2])^1.5                                            ; distance to the body center cubed for each packet (units of km)
    r_Sun3  = (sun_particle[0,*]^2 + sun_particle[1,*]^2 + sun_particle[2,*]^2)^1.5  ; distance from the Sun cubed (units of km)
    helio_dist = reform(r_Sun3^(1./3.))

  ; Calculate radiation acceleration (the radaccel.pro inputs are velocity in m/s and range to be in AU -> convert from km units)
    radaccel, Line_data.line, v*1000., helio_dist / 149597871., Line_data.wavelength, Line_data.intensity, arad ; arad is in cm/2^s
    ; Note that this result is ~15% less than Smyth models for Na

  ; Scale radiation acceleration by the fraction of sunlight that each particle sees & convert from cm/s^2 to km/s^2
    arad = float(arad) * illumination(time, x, y, z) / 1.e5 

  if keyword_set(GM_Parent) then begin

    ; Find the position of the Body WRT the Parent Body for each particle's time
      CSPICE_SPKPOS, body, ephemeris_time - Obs_Body_Ltime - REFORM(time), 'J2000', 'NONE', Parent_ID, Parent_At_Step, ltime
      Parent_particle = Parent_At_Step[0:2,*] + [x, y, z]                                         ; Vectors from the Parent Body to the particles (units of km)
      r_Parent3 = (Parent_particle[0,*]^2 + Parent_particle[1,*]^2 + Parent_particle[2,*]^2)^1.5  ; Distance from the Parent Body cubed (units of km)

    ; Cartesian gravity in km/s^2
;      gravity_x = (GM_Body * x / r_Body3) + (GM_Sun * sun_particle[0,*] / r_Sun3) + (GM_Parent * parent_particle[0,*] / r_Parent3)
;      gravity_y = (GM_Body * y / r_Body3) + (GM_Sun * sun_particle[1,*] / r_Sun3) + (GM_Parent * parent_particle[1,*] / r_Parent3)
;      gravity_z = (GM_Body * z / r_Body3) + (GM_Sun * sun_particle[2,*] / r_Sun3) + (GM_Parent * parent_particle[2,*] / r_Parent3)      
      gravity_x = (GM_Body * x / r_Body3) + (GM_Parent * parent_particle[0,*] / r_Parent3)   ; Ignore solar gravity for frames in a planetary barycenter
      gravity_y = (GM_Body * y / r_Body3) + (GM_Parent * parent_particle[1,*] / r_Parent3)
      gravity_z = (GM_Body * z / r_Body3) + (GM_Parent * parent_particle[2,*] / r_Parent3) 
  endif else begin
    ; Cartesian gravity in km/s^2
      gravity_x = (GM_Body * x / r_Body3) + (GM_Sun * sun_particle[0,*] / r_Sun3)                                  
      gravity_y = (GM_Body * y / r_Body3) + (GM_Sun * sun_particle[1,*] / r_Sun3) 
      gravity_z = (GM_Body * z / r_Body3) + (GM_Sun * sun_particle[2,*] / r_Sun3) 
  endelse

  ; Note that particle orbit's will not be correct without the Sun's / Parent bodies' gravity.
  ; Project the acceleration along a vector along the sun-atom line
    radaccel_x = arad * sun_particle[0,*] / helio_dist
    radaccel_y = arad * sun_particle[1,*] / helio_dist
    radaccel_z = arad * sun_particle[2,*] / helio_dist

  ; Consider adding a small component to radiation acceleration due to surface reflection here.
    ax = gravity_x + radaccel_x                                                                                                                                    
    ay = gravity_y + radaccel_y
    az = gravity_z + radaccel_z
    
;
;    if keyword_set(debug) then begin
;      print, 'Gravitational Acceleration in m/s^2 =', mean(sqrt((gravity_x)^2+(gravity_y)^2+(gravity_z)^2)*1.e3)  ; Gravitational acceleration in m/s^2)
;      print, 'Radiation Acceleration in cm/s^2 =', mean(sqrt((radaccel_x)^2+(radaccel_y)^2+(radaccel_z)^2)*1.e5)  ; Radiation acceleration in cm/s^2
;    endif

  return, [ax,ay,az]
end

;******************************************************************************************************************************

function RK_4_step, loc, dt, ionizelife
  COMMON Gravity, GM_Body, GM_Sun, GM_Parent
  COMMON Model_shared, Body, Ephemeris_time, Obs_Body_Ltime, Parent_ID, Seed, Directory, Particle_data, Line_data, Debug

  ; See Press et al. Numerical Recipees Eq. 17.1.3

  ;loc:   10,(number of packet) array
  ;       0 = x
  ;       1 = y
  ;       2 = z
  ;       3 = r
  ;       4 = fraction of original "packet" content
  ;       5 = vx
  ;       6 = vy
  ;       7 = vz
  ;       8 = age = how long this particle has been tracked
  ;       9 = time between particle release and image taken (seconds)

  ;x,y,z in km
  ;vx,vy,vz in km/s

  t = loc[9,*] - loc[8,*]
  if n_elements(dt) gt 1 then dt = transpose(dt)

  ; a        = f( x_body, X_parent, v_radial, t )
  ; v_radial = radial component of the atoms velocities in units of km/s

  ; RK4 Step 1
    ; The XYZ J2000 state vector with respect to the parent is:
      x1 = loc[0:2,*]     ; Position with respect to the origin (parent) at whatever this body is orbiting
      v1 = loc[5:7,*]     ; Velocity with respect to the origin (parent) at whatever this body is orbiting

    ; Get the XYZ with respect to the body itself...
      cspice_spkezr, body, ephemeris_time - OBS_BODY_LTIME - reform(t), 'J2000', 'NONE', Parent_ID, Parent_Body_State1, ltime ; Get the particles position wrt body
      X_WRT_Body1   = x1 - Parent_Body_State1[0:2,*]
      V_WRT_Body1   = v1 - Parent_Body_State1[3:5,*]

      cspice_spkezr, body, ephemeris_time - OBS_BODY_LTIME - reform(t), 'J2000', 'NONE', 'Sun', Body_Sun_state1, ltime ; Get the body's position wrt the Sun SHOULD HAVE LIGHT TIME CORRECTION
      X_WRT_Sun1 = Body_Sun_state1[0:2,*] + X_WRT_Body1  
      V_WRT_Sun1 = Body_Sun_state1[3:5,*] + V_WRT_Body1  
      
      V_radial1  = TOTAL( X_WRT_Sun1 * V_WRT_Sun1, 1 ) / sqrt(TOTAL( X_WRT_Sun1^2, 1 ) ) ;particle radial velocity in km/s
      a1 = RK_4_acceleration(X_WRT_Body1, X_WRT_Sun1, t, v_radial1)
      ;a1 = RK_4_acceleration(X_WRT_Body1, X1, t, v_radial1)
  
 ; sun_pre_step = X_WRT_Sun1 ; X,Y,Z particles positions wrt Sun before the timestep (used for photoionization)

  ; RK4 Step 2
    ; The XYZ with respect to the parent is just:
      X2 = X1 + .5*[dt,dt,dt]*v1 ; Position with respect to the origin at whatever this body is orbiting
      V2 = V1 + .5*[dt,dt,dt]*a1 ; Velocity with respect to the origin at whatever this body is orbiting

    ; Get the XYZ with respect to the body itself...
      cspice_spkezr, body, ephemeris_time - OBS_BODY_LTIME - (reform(t)- .5*dt), 'J2000', 'NONE', Parent_ID, Parent_Body_State2, ltime ; Get the particles position wrt body
      X_WRT_Body2   = x2 - Parent_Body_State2[0:2,*]
      V_WRT_Body2   = v2 - Parent_Body_State2[3:5,*]

      cspice_spkezr, body, ephemeris_time - OBS_BODY_LTIME - (reform(t)- .5*dt), 'J2000', 'NONE', 'Sun', Body_Sun_state2, ltime ; Get the body's position wrt the Sun SHOULD HAVE LIGHT TIME CORRECTION
      X_WRT_Sun2 = Body_Sun_state2[0:2,*] + X_WRT_Body2
      V_WRT_Sun2 = Body_Sun_state2[3:5,*] + V_WRT_Body2

      V2_radial  = TOTAL( X_WRT_Sun2 * V_WRT_Sun2, 1 ) / sqrt(TOTAL( X_WRT_Sun2^2, 1 ) )
      ;a2 = RK_4_acceleration(X_WRT_Body2, X2, t - .5*dt, v2_radial)
      a2 = RK_4_acceleration(X_WRT_Body2, X_WRT_Sun2, t - .5*dt, v2_radial)

  ; RK4 Step 3
    ; The XYZ with respect to the parent is just:
      X3 = X1 + .5*[dt,dt,dt]*v2 ; Position with respect to the origin at whatever this body is orbiting
      V3 = V1 + .5*[dt,dt,dt]*a2 ; Velocity with respect to the origin at whatever this body is orbiting

    ; Get the XYZ with respect to the body itself...
      ;cspice_spkezr, body, ephemeris_time - (reform(t)- .5*dt), 'J2000', 'NONE', Parent_ID, Parent_Body_State3, ltime ; Get the particles position wrt body
      X_WRT_Body3   = x3 - Parent_Body_State2[0:2,*] ;since Parent_Body_State2 = Parent_Body_State3
      V_WRT_Body3   = v3 - Parent_Body_State2[3:5,*] ;since Parent_Body_State2 = Parent_Body_State3

      ;cspice_spkezr, body, ephemeris_time - (reform(t)- .5*dt), 'J2000', 'NONE', 'Sun', Body_Sun_state3, ltime ; Get the body's position wrt the Sun SHOULD HAVE LIGHT TIME CORRECTION
      X_WRT_Sun3 = Body_Sun_state2[0:2,*] + X_WRT_Body2 ; since Body_Sun_state2 = Body_Sun_state3
      V_WRT_Sun3 = Body_Sun_state2[3:5,*] + V_WRT_Body2 ; since Body_Sun_state2 = Body_Sun_state3

      V3_radial  = TOTAL( X_WRT_Sun3 * V_WRT_Sun3, 1 ) / sqrt(TOTAL( X_WRT_Sun3^2, 1 ) )
      ;a3 = RK_4_acceleration(X_WRT_Body3, X3, t - .5*dt, v3_radial)
      a3 = RK_4_acceleration(X_WRT_Body3, X_WRT_Sun3, t - .5*dt, v3_radial)
      
  ; RK4 Step 4
    ; The XYZ with respect to the parent is just:
      X4 = X1 + [dt,dt,dt]*v3 ; Position with respect to the origin at whatever this body is orbiting
      V4 = V1 + [dt,dt,dt]*a3 ; Velocity with respect to the origin at whatever this body is orbiting
      
    ; Get the XYZ with respect to the body itself...
      cspice_spkezr, body, ephemeris_time - OBS_BODY_LTIME - (reform(t)- dt), 'J2000', 'NONE', Parent_ID, Parent_Body_State4, ltime ; Get the particles position wrt body
      X_WRT_Body4   = x4 - Parent_Body_State4[0:2,*]
      V_WRT_Body4   = v4 - Parent_Body_State4[3:5,*] 
      
      cspice_spkezr, body, ephemeris_time - OBS_BODY_LTIME - (reform(t)- dt), 'J2000', 'NONE', 'Sun', Body_Sun_state4, ltime ; Get the body's position wrt the Sun SHOULD HAVE LIGHT TIME CORRECTION
      X_WRT_Sun4 = Body_Sun_state4[0:2,*] + X_WRT_Body4 ; since Body_Sun_state2 = Body_Sun_state3
      V_WRT_Sun4 = Body_Sun_state4[3:5,*] + V_WRT_Body4 ; since Body_Sun_state2 = Body_Sun_state3

      V4_radial  = TOTAL( X_WRT_Sun4 * V_WRT_Sun4, 1 ) / sqrt(TOTAL( X_WRT_Sun4^2, 1 ) )
      ;a4 = RK_4_acceleration(X_WRT_Body4, X4, t - dt, v4_radial)
      a4 = RK_4_acceleration(X_WRT_Body4, X_WRT_Sun4, t - dt, v4_radial)

  ; And the result over this time step...

  ; Scale photo-ionization with 1 / distance^2 from the Sun and the fractional illumination
    check_illum = illumination(t, X_WRT_Body1[0,*], X_WRT_Body1[1,*], X_WRT_Body1[2,*]);check the fraction of the solar disc seen by the atoms (units of km)

  ; The e-folding photo-ionization lifetime in at this distance WRT the Sun is [units of seconds]:
    photlife   = ((sqrt(X_WRT_Sun1[0,*]^2 + X_WRT_Sun1[1,*]^2 + X_WRT_Sun1[2,*]^2) / 149597871.)^2) * (ionizelife/check_illum)

    loc[4,*]   = loc[4,*] * exp(-dt / photlife) ; Ionize all tracked atoms over the timestep
    loc[0:2,*] = loc[0:2,*] + ([dt,dt,dt]/6.)*(v1 + 2.*v2 + 2.*v3 + v4)
    loc[5:7,*] = loc[5:7,*] + ([dt,dt,dt]/6.)*(a1 + 2.*a2 + 2.*a3 + a4)
    loc[8,*]   = loc[8,*] + dt                  ; Advance by one timestep
    t          = loc[9,*] - loc[8,*]
  return, loc
end

pro RK4_integrate_adaptive, loc, reimpact_loc, bounce, ionizelife, atoms_per_packet, timestep, step_type, thermal_accom_coeff = thermal_accom_coeff
  COMMON Model_shared, Body, Ephemeris_time, Obs_Body_Ltime, Parent_ID, Seed, Directory, Particle_data, Line_data, Debug
  COMMON Output_shared, Plot_range, Output_Size_In_Pixels, Output_Title, Center_in_frame, viewpoint, FOV, N_ticks, Tickstep, Observatory, Above_Ecliptic, Boresight_Pixel, Aperture_Corners
  COMMON Gravity, GM_Body, GM_Sun, GM_Parent

  t = loc[9,*] - loc[8,*]                                                 ; time remaining to track particle

  platescale = float(float(Output_Size_In_Pixels[0])/(float(plot_range))) ; plate scale of the output in pixels per body radii

  ; Get the relevant gravitational constants, Body, Sun, and, if the Body is a moon, Parent Body

  ; Sun
    cspice_bodvrd, 'Sun', "GM", 3, GM_Sun ; The Keplerian GM of the Sun in units km^3/s^2
    GM_Sun = float(-GM_Sun[0])            ; negate (attractive)

  ; Body
    cspice_bods2c, body, bodyID_code, found
    radius_found = cspice_bodfnd( bodyID_code, "RADII" )
    Mass_found = cspice_bodfnd( bodyID_code, "GM" )
    if (radius_found and mass_found) then begin
      cspice_bodvrd, Body, 'RADII', 3, Body_radius                                                        ; Find the body's radius in Km
      flat = (Body_radius[0] - Body_radius[2])/Body_radius[0]
      Body_radius = float(Body_radius[0])
      cspice_bodvrd, body, "GM", 3, GM_Body                                                               ; Find G*Mass in in Km^3/s^2
      if keyword_set(debug) then print, 'Body Gravitational Constant: G * Mass (m^3/s^2) =', GM_Body*1.e9 ; Mercury is -2.20356e13 m^3/s^2
      GM_Body = float(-GM_Body[0])                                                                        ; Negate (attractive)
      GM_Parent = GM_Sun
    endif else begin                                                                                      ; Error Handling
      VAR = strcompress('BODY'+string(bodyID_code)+'*', /remove_all)
      cspice_gnpool, VAR, 0, 10, 81, kervar, found
      print, 'Insufficient data in the SPICE kernal pool needed to calculate gravitational constants. Available information:'
      if (found) then begin
        for i=0, n_elements(kervar)-1 do begin
          print, '   Variable ' + string(i) + ' matching ' + VAR $
            + ' : ', kervar[i]
        endfor
      endif else begin
        print, 'Failed to find  ' + VAR + ' in the kernel pool'
      endelse
      if (strmid(body, 0, 3) eq '100') then begin
        GM_Body = -1.e-6                                                                                  ; Use a small gravitational constant for small bodies
        Body_radius = 1.
        GM_Parent = GM_Sun
      endif else stop                                                                   ; Solar Sytem Barycenter is the "parent" that these bodies orbits, unless moon      
    endelse

  ; Parent body (Moons only)
    if ((STRPOS(string(bodyID_code), '99') eq -1) and (strlen(strcompress(bodyID_code, /remove_all)) eq 3)) then begin ; if it's a Moon then ....
      Parent_Planet = strcompress(strmid(strcompress(bodyID_code, /remove_all),0,1) + '99')
      cspice_bodvrd, Parent_Planet, "GM", 3, GM_Parent                                                                  ; Find G*Mass in in Km^3/s^2
      if keyword_set(debug) then print, 'Parent Body Gravitational Constant: G * Mass (m^3/s^2) =', GM_Parent*1.e9 
      GM_Parent = -GM_Parent[0]                                                                                     ; negate (attractive)
    endif

  ; Deal with surface interactions. If so we'll also need a surface temperature model, run it if it hasn't been run already
    IF KEYWORD_SET(bounce) THEN BEGIN                             
      Thermal_model_exists    = File_search(Directory+body+'_Temperature.nc')
      if Thermal_model_exists eq Directory+body+'_Temperature.nc' then begin
        ID            = NCDF_OPEN(Thermal_model_exists, /NOWRITE)
        tempID        = ncdf_varid(ID, 'temp')
        latID         = ncdf_varid(ID, 'lat')
        lonID         = ncdf_varid(ID, 'lon')
        if body eq 'Mercury' then orbital_lonID = ncdf_varid(ID, 'true_anomaly') else orbital_lonID = ncdf_varid(ID, 'orbital_lon') ; Mercury uses true anomaly angle
        ncdf_varget, ID, tempID, Temperature
        ncdf_varget, ID, latID, Thermal_lat
        ncdf_varget, ID, lonID, Thermal_lon         ; lon increases with time
        ncdf_varget, ID, orbital_lonID, orbital_lon ; orbital Longitude, or in Mercury's case, true anomaly angle (zero eq Jovian Noon)
        ncdf_close, ID
      endif else Surface_Temperature, Directory, Body, Time, Temperature_map, Silent = Silent ; Fix needed, the time input here might need to be an orbital longitude/TAA of 0.0 (unsure of this) 
      reimpact_loc = []                              ; see definitions in reimpacts array within while loop below
    endif  

  ; Tolerances for the adaptive RK step      
    ; these could be adaptive themselves i.e.
    ; tolspace = platescale^(-1.)                  ; Spatial tolerence is the size of a pixel in units of body radii
    ; tolvelo  = tolspace / t                      ; vector: the more time remaining, the smaller the v error has to be
    ; time_error = body_radius * tolspace / 100.   ; time error allowed - time for 100 km/s particle to travel tolspace
    ; Keep it simple here:    
      tolspace   = 25.   ; km   SPATIAL TOLERANCE PER STEP
      tolvel     = 0.25  ; km/s VELOCITY TOLERANCE PER STEP  
      time_error = 10.   ; s NOT ACTUALLY UTILIZED (ACCEPT FOR THE FINAL LOOP TIMEOUT)

  sl     = size(loc)
  h      = float(replicate(timestep, sl[2]))  ; Initial guess for a time step size (seconds)
  safety = .9
  shrink = -.25
  grow   = -.2
  fcor   = 1./15.                             ; Fifth-order correction factor per numerical recipees. 
  done   = 0
  count  = 0
  hold   = h
  
  max_integration = max(t)

  ;***********************************************
  while not(done) do begin

    ; Arrays are full size here
      t = loc[9,*] - loc[8,*]                                       ; Time remaining to track particle (before the RK timestep)
      moretogo = where(t gt 0., /NULL)                              ; Indicies of particles still being integrated
      if moretogo eq !Null then begin
        done = 1
        BREAK                                                       ; Break the while loop if all the integrations are completed
      endif

    ; Now only include particles still being tracked
      h = hold[moretogo]
      t = t[moretogo]
      h = (h le t)*h + (h gt t)*t                                   ; h cannot be bigger than the integration time left, clamp it there if need be

    Case Step_Type of                                               ; Is the integration timestep the same for all particles, or adaptive (per Numerical Recipes) 
      'Fixed': loc[*,moretogo] = RK_4_step(loc[*,moretogo], h, ionizelife)   
      'Adaptive': begin
        
        ; Monitor the time step distribution for these tolerances
          if keyword_set(debug) then begin 
            mm = minmax(h)
            if mm[0] ne mm[1] then cghistoplot, h/60., mininput = timestep/120., maxinput = max_integration/60., xtitle = 'Timestep (minutes) --- axis maximum is the max integration time'
          endif
        
        ; Do two half-steps
          loc_HS = loc[*,moretogo]
          loc_HS = RK_4_step(loc_HS, 0.5*h, ionizelife) 
          loc_HS = RK_4_step(loc_HS, 0.5*h, ionizelife) 

        ; Do one whole step
          loc_WS = loc[*,moretogo]
          loc_WS = RK_4_step(loc_WS, h, ionizelife) 

        ; Compare
          delta            = loc_HS - loc_WS                         ; Difference in result for two time steps (more accurate minus less accurate)
          Abs_delta        = abs(loc_HS - loc_WS)                    ; Abs_delta is now the abs. value of the difference in loc for the two time steps
          Abs_delta[0:2,*] = Abs_delta[0:2,*] / tolspace             ; The absolute differences relative to the tolerated spatial errors
          Abs_delta[5:7,*] = Abs_delta[5:7,*] / tolvel               ; The absolute differences relative to the tolerated velocity errors 
          
          errmax = reform(Abs_delta[0,*])
          errmax = (errmax ge Abs_delta[1,*])*errmax + (errmax lt Abs_delta[1,*])*Abs_delta[1,*]
          errmax = (errmax ge Abs_delta[2,*])*errmax + (errmax lt Abs_delta[2,*])*Abs_delta[2,*]
          errmax = (errmax ge Abs_delta[5,*])*errmax + (errmax lt Abs_delta[5,*])*Abs_delta[5,*]
          errmax = (errmax ge Abs_delta[6,*])*errmax + (errmax lt Abs_delta[6,*])*Abs_delta[6,*]
          errmax = (errmax ge Abs_delta[7,*])*errmax + (errmax lt Abs_delta[7,*])*Abs_delta[7,*]
          ; >1 is too much, <1 is too little
        
          loc1    = loc_HS - loc[*,moretogo]                         ; loc1 is now difference between new (2 half-step) and old values
          diffmax = reform(abs(loc1[0,*]))
          diffmax = (diffmax ge abs(loc1[1,*]))*diffmax + (diffmax lt abs(loc1[1,*]))*abs(loc1[1,*])
          diffmax = (diffmax ge abs(loc1[2,*]))*diffmax + (diffmax lt abs(loc1[2,*]))*abs(loc1[2,*])
          
          htemp   = h
          zererr  = where(errmax le 0)
          if zererr[0] ge 0 then begin
            errmax[zererr] = 1.
            h[zererr]      = htemp[zererr] * 10. ; if any maximum errors are zero or negative grow the timestep by 10 x  
          endif

        ; Adaptively change the timestep, cf. Numerical Recipees 16.2 (page 719 in the 1992 ISBN 0-521-43108-5 version)   
          h       = safety*htemp*(errmax^(grow*(errmax le 1)+shrink*(errmax gt 1)))
    
        ; Some of those timestep changes could have changed too much, clamp down and nbound the changes by 0.1 and 5.0     
          shrank_too_much   = where( (h / htemp) lt 0.1, /NULL, N_shrank_too_much )    ; no more than a factor of 10 stepsize decrease
          if N_shrank_too_much gt 0 then h[shrank_too_much] = 0.1 * htemp[shrank_too_much]
          grew_too_much     = where( (h / htemp) gt 5.0, /NULL, N_grew_too_much )      ; no more than a factor of 5 stepsize increase
          if N_grew_too_much gt 0 then  h[grew_too_much] = 5.0 * htemp[grew_too_much]

          loc_orig = loc[*,moretogo]                                ; the original values 
          loc[*,moretogo] = loc1 + loc[*,moretogo] + delta*fcor     ; *** loc now has new values using 2 halfsteps ***

        ; now for the elements where errors are too big
          toobig = where(errmax gt 1.)
          if toobig[0] ge 0. then loc[*,moretogo[toobig]] = loc_orig[*,toobig]  ; reset to original values and don't take these steps

        ; what if no change over the step happened at all ?
          insig = where(diffmax eq 0.)
          if insig[0] ge 0. then begin
            h[insig] = h[insig]*10.
            calldone = where(h[insig]-t[moretogo[insig]] gt 0, N_calldone, /NULL)       
            if N_calldone ge 0 then loc[8,moretogo[insig[calldone]]] = loc[9,moretogo[insig[calldone]]]  ; These integrations are almost done, but the remaining timestep is too short to make any difference
          endif ; no change in position
          
      end ; End Adaptive Step Case (versus fix time step case. 
    endcase

    
    t = loc[9,*] - loc[8,*]                                 ; Time remaining to track particle
    
    ; Check for particles hitting the body after this time step  
      cspice_spkpos, body, ephemeris_time - Obs_body_Ltime - reform(t[moretogo]), 'J2000', 'NONE', Parent_ID, body_coords, light_time ; body_coords is body's [x,y,z] (km) WRT Parent
      loc[3,moretogo] = sqrt((loc[0,moretogo]-body_coords[0,*])^2. + (loc[1,moretogo]-body_coords[1,*])^2. + (loc[2,moretogo]-body_coords[2,*])^2.)
      Hit_Body = where(loc[3,moretogo] lt body_radius[0], N_reimpacts)
      if Hit_Body[0] ne -1 then begin
        IF KEYWORD_SET(bounce) THEN BEGIN                   ; turn on the bouncing of particles off the surface if bounce = 1
          Number_that_bounced = 0

          ; calculate the surface temperature at the re-impact location
          ; transform landing coordinates to body-fixed cartesian and then to planetographic coordinates   
            cspice_sxform, 'J2000', 'IAU_' + body, ephemeris_time - t[Hit_Body], J2000_to_body_fixed_xform 
            
            BF_State       = Dblarr(6, N_reimpacts)         ; needed? 
            reimpacts      = Fltarr(6, N_reimpacts)         ; info on the reimpacting particles
                                                            ; 0 = solar longitude          (-180 to 180)
                                                            ; 1 = solar latitude           (-90 to 90)
                                                            ; 2 = planetographic longitude (0 to 360)
                                                            ; 3 = planetographic latitude  (-90 to -90)      
                                                            ; 4 = orbital longitude / TAA  (0 to 360)
                                                            ; 5 = fractional content       (0 to 1)
            reimpacts[5,*] = loc[4,moretogo[Hit_Body]] * atoms_per_packet                                               ; log the number of particles that reimpacted the surface over this timestep
            
            for i = 0L, long(N_reimpacts - 1) do begin
              BF_State[*,i] = transpose( J2000_to_body_fixed_xform[*,*,i] ) # loc[0:5, moretogo[Hit_Body[i]]]           ; Body-Fixed cartesian states
              cspice_recpgr, Body, BF_State[0:2,i], Body_radius[0], flat, p_lon, p_lat, p_alt                           ; Particle coords planetographic longitude definition
              
              cspice_subslr, 'Near point: ellipsoid', body, ephemeris_time - t[Hit_Body[i]], 'IAU_'+body, 'none', viewpoint, sub_solar_point_body_fixed, trgepc, srfvec 
              cspice_recpgr, Body, sub_solar_point_body_fixed, Body_radius[0], flat, subslr_lon, subslr_lat, subslr_alt ; Sub-solar point coords planetographic longitude definition
      
              ; Find Mercury's True Anomaly Angle
                cspice_spkezr, Body, ephemeris_time - t[Hit_Body[i]], 'J2000', 'LT+S', 'Sun', state, ltime
                cspice_oscelt, state, ephemeris_time - t[Hit_Body[i]], -GM_Sun, elts      ; Compute osculating orbital elements
                ecc = ELTS[1]                                                             ; Eccentricity
                MA  = ELTS[5]                                                             ; Mean Anomaly radians
                ; Now get true anomaly from Mean anomaly and Eccentricity, see https://en.wikipedia.org/wiki/True_anomaly
                ; Roy, A.E. (1988). Orbital Motion (1 ed.). Bristol, UK; Philadelphia, PA: A. Hilger. ISBN 0852743602.
                True_Anomaly = !radeg * ( MA + (2.*ecc - 0.25*ecc^3)*sin(MA) + 1.25*(ecc^2)*sin(2.*MA)+ (13./12.)*(ecc^3)*sin(3.*MA) ) ; True Anomaly Angle in degrees
                
                solar_lat = (P_lat - subslr_lat)*!radeg
                solar_lon = (P_lon - subslr_lon)*!radeg 
                if solar_lon lt -180. then solar_lon = solar_lon + 360.
                if solar_lon gt 180. then solar_lon = solar_lon - 360.

                reimpacts[0:1,i]  = [solar_lon, solar_lat]
                reimpacts[2:3,i]  = [P_lon, P_lat] * cspice_DPR()
                reimpacts[4,i]    = True_Anomaly 
              endfor

              ; Hack: I've not checked the E/W longitude convention in the thermal maps agrees with planetograph E/W. Murphy says it's wrong

              ; "Temperature" 3D thermal maps produced by surface_temperature.pro have dimensions [longitude, latitude, orbital longitude or TAA]
              ; Its number of elements in the thermal maps depends on the resolution though 
              ; Interpolate to find the fractional indices of the , then interpolate within the thermal map's 3D datacube
                solar_lambda        = interpol(findgen(n_elements(thermal_lon)), thermal_lon, reimpacts[0,*])  ; Index of re-impacted particle's longitude from the sub-solar point
                solar_mu            = interpol(findgen(n_elements(thermal_lat)), thermal_lat, reimpacts[1,*])  ; Index of re-impacted particle's latitude from the sub-solar point
                TAA_index           = interpol(findgen(n_elements(orbital_lon)), orbital_lon, reimpacts[4,*])  ; Index of the instantaneous orbital longitude / TAA at the re-impact time
                Local_surface_temps = interpolate(Temperature, solar_lambda, solar_mu, TAA_index)
             
              ; Use the local surface temperatures to find which particle bounce/stick with a Monte-Carlo technique 
                Sticking_prob       = sticking_coefficient(particle_data.name, Local_surface_temps) ; 0 to 1 probabilities for sticking, from experimental data on re-adsorption vs. temperture             
                dice                = transpose(randomu(seed, N_reimpacts))                                    ; Monte-Carlo method, roll the dice. . . 
                stick_ind           = where(dice lt Sticking_prob, N_stuck, complement = bounce_ind, Ncomplement = N_bounced, /Null) ; bounce or stick ? 

              ; For those particles that stick...
                loc[4,moretogo[Hit_Body[stick_ind]]] = 0.                                                      ; set the fractional content of packets that hit the surface to zero
                loc[8,moretogo[Hit_Body[stick_ind]]] = loc[9,moretogo[Hit_Body[stick_ind]]]                    ; set the time left in the integration to zero (because t = loc(9,*)-loc(8,*))

                ; concatenate the contribution to the surface reservoir during this time step 
                  reimpact_loc = [ [reimpact_loc], [reimpacts[*,stick_ind]] ]

              ; For those particles that bounce... 
                if N_bounced gt 0L then begin

                  ; indices within "loc" of particles within loc that we want to bounce
                    bounce_indicies = moretogo[Hit_Body[bounce_ind]]
                  
                  ; Reset the particle coordinates to just above the surface
                    loc[0:2, bounce_indicies] = loc[0:2, bounce_indicies] * 1.00001 * $
                                                body_radius[0] / rebin(loc[3, bounce_indicies], 3, N_bounced)
                                                           
                  ; Reset the particle velocities, with optional thermal acccomodation to the local surface temp     
                    V_thermal = sqrt(2.*1.3806503e-23*Local_surface_temps/particle_data.mass) * 1.e-3               ; surface thermal speed, in km/s
                    V_initial = sqrt(loc[5, bounce_indicies]^2+loc[6, bounce_indicies]^2+loc[7, bounce_indicies]^2) ; re-impact velocity
                    V_final   = thermal_accom_coeff*V_thermal+(1.-thermal_accom_coeff)*V_initial                    ; give it a new speed in km/sec WRT the surface                                    
  
                  ; Assign new trajectories, bounce trajectories are weighted a COSINE DEPENDENCE WRT the surface normal
                    r     = findgen(N_bounced) / N_bounced                        ; a zero to 1 array of n_particle increments
                    theta = fltarr(N_bounced)                                     ; the independent variable that we're weighting over, in this case the angle each particle will have
    
                    ; Define the probability distribution function that we're going to weight over.
                      distribution_function = sin(!pi*r/2.) * cos(!pi*r/2.)
    
                    sum = TOTAL( distribution_function, /CUMULATIVE )             ; get the cumulative sum of the distribution over elements 0 to i, ie, the cumulative distribution function
                    sum = sum/max(sum)                                            ; normalize the CDF to a maximum of one
    
                    for i = 1L, long(N_bounced-1) do begin                        ; step through the CDF, I don't see a way to vectorize this
                      inrange = where((r ge sum[i-1]) and (r lt sum[i]), count_ind, /null)
                      if count_ind gt 0 then theta[inrange] = i
                    endfor
                    theta = (theta/N_bounced) * (!pi/2.)                          ; scale the distribution from 0 to pi/2 radians.
                    phi   = 360.*randomu(seed, N_bounced) / !RADEG                ; bounces are isotropic in azimuth 
    
                    ; -> randomize them so that they differ for each particle following the random number seed.
                      randomNumbers = RANDOMU(seed, N_bounced)
                      scrambled_indicies = SORT(randomNumbers)
                      theta[indgen(N_bounced)] = theta[scrambled_indicies]
                     
                    loc[5,bounce_indicies] = sin(theta)*cos(phi)  ; the x-direction
                    loc[6,bounce_indicies] = sin(theta)*sin(phi)  ; the y-direction
                    loc[7,bounce_indicies] = cos(theta)           ; the z-direction
                    z_axis = [0., 0., 1.]                         ; define the z-axis in vector space.
                   
                    ; Rotate to align the z axis in the velocity vector distribution to the surface normal vectors. . .
                      for q = 0, N_bounced-1 do begin
                        Surface_Normal              = loc[0:2,bounce_indicies[q]] ; Vectors normal to the local surface.
                        Velocity_Vector_WRT_Surface = loc[5:7,bounce_indicies[q]] ; Velocity unit vectors in an arbitary Z is up frame, we want to rotation these into a Surface normal is up frame
                        cspice_vcrss, surface_normal, z_axis, surface_normal_cross_z           ; the cross product of the two. This resultant vector is the axis of rotation
                        Suface_Normal_Angle_WRT_Z = cspice_vsep(surface_normal,z_axis)         ; the angle between z axis and the surface vector where the particle is released.
  
                        ; rotate the velocity vector so that the z axis (which it was generated WRT) is now the surface normal.
                          cspice_vrotv, Velocity_vector_WRT_Surface, surface_normal_cross_z, -(Suface_Normal_Angle_WRT_Z), Velocity_vector_WRT_Absolute
                          loc[5:7,bounce_indicies[q]] = Velocity_vector_WRT_Absolute * V_final[q]   ; write the vector back into the loc array, and scale it velocity                  
                      endfor
               endif ; N_bounced gt 0
               print, string(N_reimpacts, format = '(I6)'),' particles collided with '+body+'''s surface,', $
                      string(N_bounced, format = '(I6)'),' bounced (T-depen Sticking & Thermal Accommodation = '+strcompress(thermal_accom_coeff)+')'
        ENDIF ELSE BEGIN                                              ; bounce keyword is not set
          loc[4,moretogo[Hit_Body]] = 0.                              ; set the fractional content of packets that hit the surface to zero
          loc[8,moretogo[Hit_Body]] = loc[9,moretogo[Hit_Body]]       ; set the time left in the integration to zero (because t = loc(9,*)-loc(8,*))
          print, string(N_reimpacts, format = '(I6)'),' particles collided with '+body+'''s surface (bounce keyword not set at input)'
        ENDELSE    
      endif

      hold[moretogo] = h
      ;t = loc[9,*] - loc[8,*]                                         ; time remaining to track particle old and had a bug, moved upwards in the timeline...remove? 
      count = count + 1
      if max(t) le time_error or count gt 5000 then done = 1
      
  endwhile
  return
end