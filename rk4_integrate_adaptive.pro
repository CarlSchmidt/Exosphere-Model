function RK_4_acceleration, Body_xyz, Sun_xyz, time, v

  ; Body_xyz is a 3 x N_particles vector of the J2000 particle locations with respect to Body center at 'time'
  ; Sun_xyz is a 3 x N_particles vector of the J2000 particle locations with respect to Sun center at 'time'

  COMMON Gravity, GM_Body, GM_Sun, GM_Parent
  COMMON Model_shared, Body, Ephemeris_time, Seed, Directory, Particle_data, Line_data, Debug

  ; Find the position of the body WRT the Sun for each particle's time
  ;CSPICE_SPKPOS, body, ephemeris_time - REFORM(time), 'J2000', 'None', 'Sun', sun_at_step, ltime ; Slight hack, this should have a light time correction, but that's slow.
  ;sun_at_step = float(sun_at_step)
  ;sun_particle = sun_at_step + [x, y, z]  ;Calculate the vectors from the sun to the particles
  x = Body_xyz[0,*]
  y = Body_xyz[1,*]
  z = Body_xyz[2,*]
  sun_particle = Sun_xyz

  ; Get the distance cubed from the Body and the Sun
  r_Body3 = ([x^2] + [y^2] + [z^2])^1.5                                            ; distance to the body center cubed for each packet (units of km)
  r_Sun3  = (sun_particle[0,*]^2 + sun_particle[1,*]^2 + sun_particle[2,*]^2)^1.5  ; distance from the Sun cubed (units of km)
  helio_dist = r_Sun3^(1./3.)

  ; Calculate radiation acceleration (the radaccel.pro inputs are velocity in m/s and range to be in AU -> convert from km units)
  radaccel, Line_data.line, v*1000., helio_dist / 149597871., Line_data.wavelength, Line_data.intensity, arad ; arad is in cm/2^s

  ; Check the fraction of the solar disc seen by the atoms (units of km)
  arad = float(arad) * illumination(time, x, y, z) / 1.e5 ; scale radiation acceleration by the fraction of sunlight that each particle sees
  ; convert radiation acceleration from cm/s^2 to km/s^2
  if keyword_set(GM_Parent) then begin
    cspice_bods2c, body, bodyID_code, found
    Parent_ID = strcompress(strmid(strcompress(bodyID_code, /remove_all),0,1) + '99')

    ; Find the position of the body WRT the Parent for each particle's time
    CSPICE_SPKPOS, body, ephemeris_time - REFORM(time), 'J2000', 'NONE', Parent_ID, Parent_At_Step, ltime
    Parent_particle = Parent_At_Step[0:2,*] + [x, y, z]  ;Calculate the vectors from the parent to the particles
    r_Parent3 = (Parent_particle[0,*]^2 + Parent_particle[1,*]^2 + Parent_particle[2,*]^2)^1.5  ;distance from the SSB cubed (units of km)
    ; Cartesian gravity in km/s^2
    gravity_x = (GM_Body * x / r_Body3) + (GM_Sun * sun_particle[0,*] / r_Sun3) + (GM_Parent * parent_particle[0,*] / r_Parent3)
    gravity_y = (GM_Body * y / r_Body3) + (GM_Sun * sun_particle[1,*] / r_Sun3) + (GM_Parent * parent_particle[1,*] / r_Parent3)
    gravity_z = (GM_Body * z / r_Body3) + (GM_Sun * sun_particle[2,*] / r_Sun3) + (GM_Parent * parent_particle[2,*] / r_Parent3)
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

  ;    if keyword_set(debug) then begin
  ;      gravity  = sqrt((gravity_x)^2+(gravity_y)^2+(gravity_z)^2)*1.e3     ; Gravitational acceleration in m/s^2
  ;      radaccel = sqrt((radaccel_x)^2+(radaccel_y)^2+(radaccel_z)^2)*1.e5  ; Radiation acceleration in cm/s^2
  ;      print, 'Gravitational Acceleration in m/s^2) =', mean(gravity)
  ;      print, 'Radiation Acceleration in cm/s^2) =', mean(radaccel)
  ;    endif

  return, [ax,ay,az]
end

;******************************************************************************************************************************

function RK_4_step, loc, dt, ionizelife
  COMMON Gravity, GM_Body, GM_Sun, GM_Parent
  COMMON Model_shared, Body, Ephemeris_time, Seed, Directory, Particle_data, Line_data, Debug

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

  ; a        = f( x, v_radial, t )
  ; v_radial = radial component of the atoms velocities in units of km/s

  ; RK4 Step 1
    x1 = loc[0:2,*] ;position with respect to the body
    v1 = loc[5:7,*] ;velocity with respect to the body 
    cspice_spkezr, body, ephemeris_time - reform(t), 'J2000', 'NONE', 'Sun', Body_state, ltime ; Get the particles' heliocentric velocity
    X_WRT_Sun = Body_state[0:2,*] + x1
    V_WRT_Sun = Body_state[3:5,*] + v1
    V_radial  = TOTAL( X_WRT_Sun * V_WRT_Sun, 1 ) / sqrt(TOTAL( X_WRT_Sun^2, 1 ) )
    a1 = RK_4_acceleration(x1, X_WRT_Sun, t, v_radial)

  sun_pre_step = X_WRT_Sun ; X,Y,Z particles positions wrt Sun before the timestep (used for photoionization)

  ; RK4 Step 2
    x2 = loc[0:2,*] + .5*[dt,dt,dt]*v1
    v2 = loc[5:7,*] + .5*[dt,dt,dt]*a1
    cspice_spkezr, body, ephemeris_time - reform(t - .5*dt), 'J2000', 'NONE', 'Sun', Body_state, ltime ; Get the particles' heliocentric velocity
    X_WRT_Sun = Body_state[0:2,*] + x2
    V_WRT_Sun = Body_state[3:5,*] + v2
    V_radial  = TOTAL( X_WRT_Sun * V_WRT_Sun, 1 ) / sqrt(TOTAL( X_WRT_Sun^2, 1 ) )
    a2 = RK_4_acceleration(x2, X_WRT_Sun, t - .5*dt, v_radial)

  ; RK4 Step 3
    x3 = loc[0:2,*] + .5*[dt,dt,dt]*v2
    v3 = loc[5:7,*] + .5*[dt,dt,dt]*a2
    X_WRT_Sun = Body_state[0:2,*] + x3 ; already have the body's state vector at time t 
    V_WRT_Sun = Body_state[3:5,*] + v3
    V_radial  = TOTAL( X_WRT_Sun * V_WRT_Sun, 1 ) / sqrt(TOTAL( X_WRT_Sun^2, 1 ) )
    a3 = RK_4_acceleration(x3, X_WRT_Sun, t - .5*dt, v_radial)

  ; RK4 Step 4
    x4 = loc[0:2,*] + [dt,dt,dt]*v3
    v4 = loc[5:7,*] + [dt,dt,dt]*a3 
    cspice_spkezr, body, ephemeris_time - reform(t - dt), 'J2000', 'NONE', 'Sun', Body_state, ltime ; Get the particles' heliocentric velocity
    X_WRT_Sun = Body_state[0:2,*] + x4
    V_WRT_Sun = Body_state[3:5,*] + v4
    V_radial  = TOTAL( X_WRT_Sun * V_WRT_Sun, 1 ) / sqrt(TOTAL( X_WRT_Sun^2, 1 ) )
    a4 = RK_4_acceleration(x4, X_WRT_Sun, t - dt, v_radial) ;RK_4_acceleration(x4[0,*], x4[1,*], x4[2,*], t - dt, v_radial)

  ; And the result over this time step...

  ; Scale photo-ionization with 1 / distance^2 from the Sun and the fractional illumination
    check_illum = illumination(t, loc[0,*], loc[1,*], loc[2,*]);check the fraction of the solar disc seen by the atoms (units of km)
  ;CSPICE_SPKPOS, body, ephemeris_time - REFORM(t), 'J2000', 'None', 'Sun', sun_at_step, ltime ; Slight hack, this should have a light time correction, but that's slow.
    sun_particle = sun_pre_step  ; Calculate the vectors from the sun to the particles
  ; The e-folding photo-ionization lifetime in seconds is:
    photlife   = ((sqrt(sun_particle[0,*]^2 + sun_particle[1,*]^2 + sun_particle[2,*]^2) / 149597871.)^2) * (ionizelife/check_illum)

    loc[4,*]   = loc[4,*] * exp(-dt / photlife) ; Ionize all tracked atoms over the timestep
    loc[0:2,*] = loc[0:2,*] + ([dt,dt,dt]/6.)*(v1 + 2.*v2 + 2.*v3 + v4)
    loc[5:7,*] = loc[5:7,*] + ([dt,dt,dt]/6.)*(a1 + 2.*a2 + 2.*a3 + a4)
    loc[8,*]   = loc[8,*] + dt                  ; Advance by one timestep
    t          = loc[9,*] - loc[8,*]
  return, loc
end

pro RK4_integrate_adaptive, loc, bounce, ionizelife, atoms_per_packet, timestep
  COMMON Model_shared, Body, Ephemeris_time, Seed, Directory, Particle_data, Line_data, Debug
  COMMON Output_shared, Plot_range, Output_Size_In_Pixels, Output_Title, Center_in_frame
  COMMON Gravity, GM_Body, GM_Sun, GM_Parent

  t = loc[9,*] - loc[8,*] ;time remaining to track particle

  platescale = float(float(Output_Size_In_Pixels[0])/(float(plot_range))) ;plate scale of the output in pixels per body radii

  ; Get the relevant gravitational constants, Body, Sun, and, if Body is a moon, Parent Body
  ; Note that gravitational effects of any moons on their parent bodies are not included as written, but the converse is true

  ; Body
    cspice_bods2c, body, bodyID_code, found
    radius_found = cspice_bodfnd( bodyID_code, "RADII" )
    Mass_found = cspice_bodfnd( bodyID_code, "GM" )
    if (radius_found and mass_found) then begin
      cspice_bodvrd, Body,'RADII', 3, Body_radius ;Find the body's radius in Km
      Body_radius = float(Body_radius[0])
      cspice_bodvrd, body, "GM", 3, GM_Body ;Find G*Mass in in Km^3/s^2
      if keyword_set(debug) then print, 'Body Gravitational Constant: G * Mass (m^3/s^2) =', GM_Body*1.e9 ; Mercury is -2.20356e13 m^3/s^2
      GM_Body = float(-GM_Body[0]) ;negate (attractive)
    endif else begin ; Error Handling
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
      stop
    endelse

  ; Sun
    cspice_bodvrd, 'Sun', "GM", 3, GM_Sun ; Find G*Mass in in Km^3/s^2
    GM_Sun = float(-GM_Sun[0])            ; negate (attractive)

  ; Parent body (Moons only)
    if ((STRPOS(string(bodyID_code), '99') eq -1) and (strlen(strcompress(bodyID_code, /remove_all)) eq 3)) then begin
      Parent_ID = strcompress(strmid(strcompress(bodyID_code, /remove_all),0,1) + '99')
      cspice_bodvrd, Parent_ID, "GM", 3, GM_Parent ;Find G*Mass in in Km^3/s^2
      if keyword_set(debug) then print, 'Parent Body Gravitational Constant: G * Mass (m^3/s^2) =', GM_Parent*1.e9 ;Saturn is ? m^3/s^2
      GM_Parent = -GM_Parent[0]  ;negate (attractive)
    endif

  ; Tolerances for the adaptive RK step                 
    ;tolspace = platescale^(-1.)                  ; Spatial tolerence is the size of a pixel in units of body radii
    ;tolvelo  = tolspace / t                      ; vector: the more time remaining, the smaller the v error has to be
    ;time_error = body_radius * tolspace / 100.   ; time error allowed - time for 100 km/s particle to travel tolspace
    
    tolspace   = 25.   ; km   SPATIAL TOLERANCE PER STEP
    tolvel     = 0.25  ; km/s VELOCITY TOLERANCE PER STEP  
    time_error = 10.   ; s NOT ACTUALLY UTILIZED (ACCEPT FOR THE FINAL LOOP TIMEOUT)

  sl         = size(loc)
  h          = float(replicate(timestep, sl[2]))  ; Initial guess for a time step size (seconds)
  safety     = .9
  shrink     = -.25
  grow       = -.2
  fcor       = 1./15.   ;fifth-order correction factor
  done  = 0
  count = 0
  hold  = h

  ;Step_Type = 'Fixed'
  Step_Type = 'Adaptive'

  max_integration = max(t)

  ;***********************************************
  while not(done) do begin

    ; Arrays are full size here
      t = loc[9,*] - loc[8,*]                                       ; Time remaining to track particle
      moretogo = where(t gt 0., /NULL)                              ; Indicies of particles still being integrated
      if moretogo eq !Null then moretogo[0] = 0

    ; Now only include particles still being tracked
      h = hold[moretogo]
      t = t[moretogo]
      h = (h le t)*h + (h gt t)*t                                   ; h cannot be bigger than the integration time left, clamp it there if need be

    Case Step_Type of
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
          
          ;new way
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
          endif ;no change in position
          
      end ; Adaptive Step
    endcase

    ; Check for particles hitting the body
      loc[3,moretogo] = sqrt(loc[0,moretogo]^2 + loc[1,moretogo]^2 + loc[2,moretogo]^2)
      Hit_Body = where(loc[3,moretogo] lt body_radius)
      if Hit_Body[0] ne -1 then begin
        loc[4,moretogo[Hit_Body]] = 0.                              ; set the fractional content of packets that hit the surface to zero
        loc[8,moretogo[Hit_Body]] = loc[9,moretogo[Hit_Body]]       ; set the time left in the integration to zero (because t = loc(9,*)-loc(8,*))
        print,n_elements(Hit_Body),'  Particles collided with the surface (no bouncing)'
      endif

      hold[moretogo] = h
      t = loc[9,*] - loc[8,*]                                       ; time remaining to track particle
      count = count + 1
      if max(t) le time_error or count gt 5000 then done = 1
  endwhile
  return
end