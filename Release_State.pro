pro Release_State, loc, speed_distribution, surface_distribution, speed

COMMON Model_shared, Body, Ephemeris_time, Obs_Body_Ltime, Parent_ID, Seed, Directory, Particle_data, Line_data, Debug
COMMON Output_shared, Plot_range, Output_Size_In_Pixels, Output_Title, Center_in_frame, viewpoint, FOV, N_ticks, Tickstep, Observatory, Above_Ecliptic, Boresight_Pixel, Aperture_Corners

  N_particles = N_elements(loc[0,*])
  ; Calculate the flatness coefficient and planetary radius
    
    if strmid(body, 0, 3) eq '100' then begin
      re   = 1.                               ; set cometary radii to 1 km
      flat = 0.                               ; spherical      
    endif else begin
      cspice_bodvrd, body, 'RADII', 3, radii
      re   = radii[0]
      rp   = radii[2]
      flat = ( re - rp ) / re
    endelse

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
    isotropic_release = 1 ;ejection angles randomly distributed over 2 pi steradians
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
    Temperature = 3500. ;
    PRINT, 'hack: NEED TO hard code this TEMPERATURE as an top level input!'
    STOP
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
  if STRPOS(speed_distribution, 'Sputtering_ST') eq 0 then begin
    Distribution = 'Sputtering_ST'
    SBE = 7.9 ; eV, yield = 0.000412 Na atoms / proton *** ALBITE *** sodium bearing plagioclase (See Morrissey et al. (2022) Table 2)
    distribution_parameters = {Sputtering_ST_structure, Distribution:Distribution, SBE:SBE}
    speed = generate_velocity_distribution(N_particles, distribution_parameters)
    isotropic_release = 1 ;ejection angles randomly distributed over 2 pi steradians NO ONE SEEMS TO KNOW IF THIS IS CORRECT, OR IF IT'S COS(NORMAL ANGLE) AS IN MB "FLUX"
  endif
  
  if STRPOS(speed_distribution, 'Johnson_2002') eq 0 then begin ;cf. Johnson et al. 2002 doi.org/10.1006/icar.2001.6763
    case Particle_data.name of
      'Na': begin
              U_param = 0.052
              x_param = 0.7
            end
      'K' : begin
              U_param = 0.02
              x_param = 0.25
            end
      else: begin
            print, 'Distribution undefined in Johnson et al. (2002) for species: ', Particle_data.name 
            stop
            end
    endcase 
    Distribution = 'Johnson_2002' 
    distribution_parameters = {Johnson_2002_Structure, Distribution:Distribution, U_param:U_param, x_param:x_param}
    speed = generate_velocity_distribution(N_particles, distribution_parameters)
    ; NOTE: Assumed preferentially upward as the cosine of the surface normal (as in a flux distribution), Leblanc et al. 2002 used 2cos(theta)
  endif


  if STRPOS(speed_distribution, 'Shematovich') eq 0 then begin
    Distribution = 'Shematovich'
    readcol, Directory+'Shematovich_2013_Fig5.txt', F='A,A', v1, v2, DELIMITER = ',' ;'Shematovich (2013) Figure 5'
    energy = float(v1) ;eV
    PDF = float(v2) ;cm-2, s-1, eV-1

    ; convert that to cm-2, s-1, (m/s)-1
      Velocity = sqrt(2.*energy*1.60218e-19/Particle_data.mass) ;convert hydrogen energy in eV to joules to velocity in m/s
      eV_over_m_per_s = 1. / sqrt(2.*1.*1.60218e-19/Particle_data.mass)
      PDF = PDF*eV_over_m_per_s ;convert hydrogen energy in eV to joules and then to velocity in m/s

    Temperature = 800. ;meaningless here, but needed to set a somewhat arbitrary upper bound on the randomn velocities within the CDF
    distribution_parameters = {Shematovich_Structure, Distribution:Distribution, Temperature:Temperature, PDF:PDF, Velocity:Velocity}
    speed = generate_velocity_distribution(n_elements(loc[9,*]), distribution_parameters)
    isotropic_release = 1
  endif
  
  ; Identify the sub-observer and sub-solar coordinates in longitude and latitude AT THE EPHEMERIS_TIME being simulated, account for the observer's light time
    cspice_subpnt, 'Near point: ellipsoid', body, ephemeris_time, 'IAU_'+body, 'LT', viewpoint, sub_observer_point_planet_frame, trgepc, srfvec
    cspice_recpgr, body, sub_observer_point_planet_frame, re, flat, subobserver_lon_ET, subobserver_lat_ET, subobserver_radius_ET
    
    ; Get an x,y,z subsolar point (in km from planet center) in the body-fixed frame, account for the observer's light time
      cspice_subslr, 'Near point: ellipsoid', body, ephemeris_time, 'IAU_'+body, 'LT', viewpoint, sub_solar_point_planet_frame, trgepc, srfvec

    ; Convert this sub-solar point to planetographic latitude and longitude, in the 'IAU_' body fixed frame
      cspice_recpgr, body, sub_solar_point_planet_frame, re, flat, subsolar_lon_ET, subsolar_lat_ET, subsolar_point_radius_ET
    
  
; ============================================== LAUNCH COORDINATES =================================================================================================
; Get the particle launch locations in body-centered J2000 coordinates
  CASE 1 OF
    STRMATCH(surface_distribution, 'Global', /FOLD_CASE): BEGIN               ; Generate isotropic latitudes and longitudes for all particles
        lat = asin(1.d - 2.d*randomu(seed, N_particles, /double))             ; gives a random number, -90 to 90 degrees (in radians) weighted towards the equator (more surface area per unit latitude)
        lon = 360.d*randomu(seed, N_particles, /double) /!RADEG               ; gives random longitudes 0 to 360 (in radians)
        if strmid(body, 0, 3) eq '100' then $                                 ; Small bodies can't do a complext lat & lon to cartesian coordinate transform, but planets can...
        release_points_body_fixed = transpose([[cos(lat)*cos(lon)], [cos(lat)*sin(lon)], [sin(lat)]]) else $
        cspice_pgrrec, body, lon, lat, replicate(0.D, N_particles), re, flat, release_points_body_fixed; convert these to cartesian, body fixed coords units of planetary radii here, not km
        release_points_J2000 = release_points_body_fixed / re                 ; since these are random and isotropic in x, y, z, the reference frame does not matter, use body radii coordinates
    END  
    STRMATCH(surface_distribution, 'Dayside', /FOLD_CASE): BEGIN
      release_points_J2000      = dblarr(3, N_Particles)
      release_points_Body_fixed = dblarr(3, N_Particles)

      ss_lon = fltarr(N_Particles) & ss_lat = fltarr(N_Particles) 
      for i = 0, N_particles - 1 do begin ; sub-solar point moves depending on the release time, so we need a for loop to do this... 
        
        ; Get an x,y,z subsolar point (in km from planet center) in the body-fixed frame at the time of release, account for the observer's light time
          cspice_subslr, 'Near point: ellipsoid', body, ephemeris_time - loc[9,i], 'IAU_'+body, 'LT', viewpoint, sub_solar_point_planet_frame, trgepc, srfvec 

        ; Convert this sub-solar point to planetographic latitude and longitude, in the 'IAU_' body fixed frame
          cspice_recpgr, body, sub_solar_point_planet_frame, re, flat, sub_solar_lon, sub_solar_lat, sub_solar_point_radius        
          
        ; Generate a randomn latitude and longitude so that the zenith of the hemisphere is along the positive y direction.
          lat = asin(1.0-2.0*randomu(seed))                            ; gives a random number between -90 and 90 degrees (in radians) 
                                                                       ; weighted towards the equator (more surface area per unit latitude)
          lon = ((180.*randomu(seed)))/!radeg - !pi/2. + sub_solar_lon ; gives a random longitude -90 and 90 degrees from the sub-solar point(in radians)

          ss_lat[i] = lat                                              ; We may want to have the packets' sub-solar coordinates...
          ss_lon[i] = Lon - sub_solar_lon                              ; sub-solar coordinates are used to weight by the solar zenith angle

        ; Convert the release point from planetographic coordinates to rectangular coordinates, still in the body fixed frame
          cspice_pgrrec, body, lon, lat, 0., re, flat, release_point_body_fixed

        ; get the rotation matrix to J2000 frame the rest of the model runs in.
          cspice_pxform, 'IAU_'+body, 'J2000', ephemeris_time - Obs_Body_Ltime - loc[9,i], body_fixed_to_J2000_rotation_matrix
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
    STRMATCH(surface_distribution, 'Point*', /FOLD_CASE): BEGIN
       print, 'I believe the light time correction needs to be applied here... do this! (Obs_Body_Ltime)'
      ;stop
        release_points_J2000      = dblarr(3, N_Particles)
        release_points_Body_fixed = dblarr(3, N_Particles) 
         
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
          release_points_J2000[*,i]      = release_point_J2000_frame / re  ; the initial packet starting points in cartesian units of planetary radii, this get's scaled back into km later
          release_points_body_fixed[*,i] = release_point_body_fixed
        endfor

    END
    STRMATCH(surface_distribution, 'From_Map*', /FOLD_CASE): BEGIN
      
      ; The technique here is to simulate random spatially isotropic points
      ; A map is then used to weight the loc[4,*] parameter which controls how many atoms are in each test particle
      
      ; Use Paul Szabo's plasma precipitation maps of Mercury during an ICME ...
      ICME = 1
      if ICME then begin
        map_filename = directory+'\Surface_reservoir\H+_Precipitation_Map_400nPa_ICME.txt'
        Proton_precip = READ_ASCII( map_filename, COMMENT_SYMBOL='#' )  
        
        map_filename = directory+'\Surface_reservoir\He++_Precipitation_Map_400nPa_ICME.txt'
        Alpha_precip = READ_ASCII( map_filename, COMMENT_SYMBOL='#' )
        
        map        = rebin(Proton_precip.field01, 360., 180.)                        ; rebin an arbitrary sized map into degree bins in lat, lon. 
        map_alphas = rebin(Alpha_precip.field01, 360., 180.)                         ; rebin an arbitrary sized map into degree bins in lat, lon. 
        
        ;Szabo_latitude_Bins = [-87.5,-82.5,-77.5,-72.5,-67.5,-62.5,-57.5,-52.5,-47.5,-42.5,-37.5,-32.5, $
        ;                       -27.5,-22.5,-17.5,-12.5, -7.5, -2.5,  2.5,  7.5, 12.5, 17.5, 22.5, 27.5, $
        ;                        32.5, 37.5, 42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5, 77.5, 82.5, 87.5]
        ;Szabo_longitude_Bins = [-177.5,-172.5,-167.5,-162.5,-157.5,-152.5,-147.5,-142.5,-137.5,-132.5, $
        ;                        -127.5,-122.5,-117.5,-112.5,-107.5,-102.5, -97.5, -92.5, -87.5, -82.5, $
        ;                         -77.5, -72.5, -67.5, -62.5, -57.5, -52.5, -47.5, -42.5, -37.5, -32.5, $
        ;                         -27.5, -22.5, -17.5, -12.5,  -7.5,  -2.5,   2.5,   7.5,  12.5,  17.5, $
        ;                          22.5,  27.5,  32.5,  37.5,  42.5,  47.5,  52.5,  57.5,  62.5,  67.5, $
        ;                          72.5,  77.5,  82.5,  87.5,  92.5,  97.5, 102.5, 107.5, 112.5, 117.5, $
        ;                         122.5, 127.5, 132.5, 137.5, 142.5, 147.5, 152.5, 157.5, 162.5, 167.5, $
        ;                         172.5, 177.5]         
        
        if keyword_set(debug) then begin
          
          ; Calculate the global precipitation rate from a map of precipitation per unit area with degree binning in lat and lon:
            latitude_steps = (findgen(181) - 90.)/!radeg
            surface_area_per_lat_lon_bin = fltarr(360, 180)
            for i = 0, 180 - 1 do begin
              surface_area_per_lat_lon_bin[*,i] = (!pi / 180.) * (re*1.e5)^2 * (sin(latitude_steps[i+1]) - sin(latitude_steps[i])) ; Area = pi/180 * R² * (sin φ1 − sin φ2) (θ1 − θ2).
            endfor

            Print, 'Global Proton Precipitation per Second =', total(surface_area_per_lat_lon_bin * map)
            Print, 'Global Alpha  Precipitation per Second =', total(surface_area_per_lat_lon_bin * map_alphas)

            ;-------------------------------- YIELDS --------------------------------------- 
               Proton_Yield = 0.000412          ; Albite value from Morrissey et al. 2022, Table 2, needs to be consistent with the SBE (7.9eV) listed in the beginning of this program. 
               Alpha_Yield  = Proton_Yield * 8. ; Vorburger et al. (2014), Wurz et al. (2007) 
            ;------------------------------------------------------------------------------ 
            
            Upward_flux_at_exobase_via_Szabo = total(surface_area_per_lat_lon_bin * map) * Proton_Yield + $
                                               total(surface_area_per_lat_lon_bin * map_alphas) * Alpha_Yield 
            Print, 'Global Sputtered Atoms/s (number to place at simulation input) =', Upward_flux_at_exobase_via_Szabo

          cgPS_Open, filename = directory+'\Surface_reservoir\Szabo_protons.eps', /ENCAPSULATED, xsize = 7.5, ysize = 7.5
          !P.font = 1
          device, SET_FONT = 'Helvetica Bold', /TT_FONT
          !p.charsize = 1.5
            loadgood
            ;cgLoadCT, 27
            axis_format = {xtickinterval:30, ytickinterval:30, ytitle:'Subsolar longitude [degrees]', xtitle:'Subsolar latitude [degrees]'}
            cgimage, map, /keep_aspect, /axes, xr = [-180, 180], yr = [-90,90], axkey=axis_format, POSITION=[0.15, 0.50, 0.75, 0.95], minval = min(map), maxval = max(map)
            cgCOLORBAR, /right, /vertical, title = 'Proton Precipition [cm!U-2!N s!U-1!N]', range = minmax(Proton_precip.field01), /fit 
            cgimage, map_Alphas, /keep_aspect, /axes, xr = [-180, 180], yr = [-90,90], axkey=axis_format, POSITION=[0.15, 0.05, 0.75, 0.55], /noerase, minval = min(map_Alphas), maxval = max(map_Alphas)
            cgCOLORBAR, /right, /vertical, title = 'Alpha Precipition [cm!U-2!N s!U-1!N]', range = minmax(Alpha_precip.field01), /fit
          cgPS_Close
        endif
        
      endif else begin
        
        ; Use the Reimpacting flux fits file ...      
          map_filename = directory+'Surface_reservoir\Mean_Reimpacting_Flux_Map.fit'
          map = mrdfits(map_filename) ; TAKE CARE THAT LONGITUDE IS IN THE PROPER CONVENTION HERE: planetographic WEST longitude is the x axis (left-handed longitude)      
      endelse

      ; Get an x,y,z subsolar point (in km from planet center) in the body-fixed frame at the time of release.
        sub_solar_point_planet_frame = fltarr(3, N_Particles)
        for i = 0, N_Particles - 1 do begin
          cspice_subslr, 'Near point: ellipsoid', body, ephemeris_time - loc[9,i], 'IAU_'+body, 'none', viewpoint, sub_solar_point_BF, trgepc, srfvec
          sub_solar_point_planet_frame[*,i] = sub_solar_point_BF
        endfor
      
      ; Convert this sub-solar point to planetographic latitude and longitude, in the 'IAU_' body fixed frame
        cspice_recpgr, body, sub_solar_point_planet_frame, re, flat, sub_solar_lon, sub_solar_lat, sub_solar_point_radius
        
      ; Generate isotropic release points in the body-fixed frame
        lat       = asin(1.d - 2.d*randomu(seed, N_particles, /double))   ; gives a random number, -90 to 90 degrees (in radians) weighted towards the equator (more surface area per unit latitude)
        lon       = 360.d*randomu(seed, N_particles, /double) /!RADEG     ; gives random longitudes 0 to 360 (in radians)

      ; Find all their sub-solar lats and lons at the release time  
        ss_lat    = lat - sub_solar_lat                                   ; We need to have the packets' sub-solar coordinates...
        ss_lon    = Lon - sub_solar_lon                                   ; Sub-solar coordinates are used to locate position relative to the precipitation map that uses MSO coords
      
      ; Convert coords to cartesian, body-fixed coords units of planetary radii here, not km
        cspice_pgrrec, body, lon, lat, replicate(0.D, N_particles), re, flat, release_points_body_fixed

      ; Convert planetographic coodinates to the J2000 frame at the ephemeris time for each particles release
        release_points_J2000      = dblarr(3, N_Particles)
        for i = 0, N_Particles - 1 do begin
          cspice_pxform, 'IAU_'+body, 'J2000', ephemeris_time - loc[9,i], frame_rotate_matrix
          cspice_mxv, frame_rotate_matrix, release_points_body_fixed[*,i], release_point_J2000_frame
          release_points_J2000[*,i]      = release_point_J2000_frame / re     ; the initial packet starting points in cartesian units of planetary radii
        endfor
        
       
     ; ------------------------------------------------------------------------------ READ ME ----------------------------------------------------------------
     ;
     ; The approach here is to generate a globally isotropic release rate 
     ; Then, the loc[4,*] fractional packet content is used to set weights above or below unity
     ; The average weight MUST remain unity to conserve all calculated quantities from the "upward flux at exobase" variable
     ; Weights are set according to the packet's initial lat & lon location relative to a unity normalized precipitation map.  
     ; Note that since packets are distributed isotropically, there are already more packets at lower latitudes, so we needed not weight for surface area (else we'd be doing that twice)  
     ; 
     ;--------------------------------------------------------------------------------------------------------------------------------------------------------

       Normalized_precipitation_map = map / mean(map)       
       loc[4,*]                     = interpolate(normalized_precipitation_map, (ss_lon*!radeg + 180.) mod 360., ss_lat*!radeg + 90.)  ; weighted fractional content of the packets
       scale_N_Packets_conserved    = N_particles / total((loc[4,*])) 
       loc[4,*]                     = loc[4,*] * scale_N_Packets_conserved                                     ; force weighted fractional content of the packets to an average of 1. 

    END  
    ELSE: PRINT, 'Undefined spatial distribution requested as an input'  
  ENDCASE

  ; THE FINAL LAUNCH COORDINATES. UNITS OF BODY-RADII. BODY-CENTERED J2000 COORDINATES
    loc[0:2,*] = temporary(release_points_J2000) ; Allocate the release coordinates to the initialize the loc array. 
                                                 ; Trajectory integrations will be done in J2000
                                                 ; These coords are later converted from planetary radii to km in the main level generic_model.pro

  if keyword_set(debug) then begin
    angles_from_normal = fltarr(N_particles) ; for inspecting the ejection angle wrt local surface normal vector

    if strmid(body, 0, 3) eq '100' then cspice_reclat, release_points_body_fixed, radius, longitudes, latitudes else $
      cspice_recpgr, body, release_points_body_fixed, re, flat, longitudes, latitudes, radius ; Convert to planetographic latitude and longitude, in the 'IAU_' body fixed frame
    window, 1
    cgplot, (longitudes*!radeg + 360.) mod 360., latitudes*!radeg, psym=3, xr=[0,360], yr=[-90,90], ytitle = 'Latitude', $             ; particle launch coordinates 
        title='Particle Launch Coordinates', xtickinterval = 30., ytickinterval = 30., xtitle = 'Planetographic W Longitude'
    cgplot, (subsolar_lon_ET*!radeg + 360.) mod 360., subsolar_lat_ET*!radeg, psym = 16, symsize = 2, color = 'red', /overplot         ; the sub-solar point, at simulation end time (not release times)
    cgplot, (subobserver_lon_ET*!radeg + 360.) mod 360., subobserver_lat_ET*!radeg, psym = 16, symsize = 2, color = 'blue', /overplot  ; the sub-observer point, at simulation end time (not release times)
    
    cgtext, (subobserver_lon_ET*!radeg + 360.) mod 360., subobserver_lat_ET*!radeg+10., 'Sub-'+viewpoint+' Point', color = 'blue', align = 0.5, charthick = 2
    cgtext, (subsolar_lon_ET*!radeg + 360.) mod 360., subsolar_lat_ET*!radeg+10., 'Sub-Solar Point', color = 'red', align = 0.5, charthick = 2
    
    
    ; Make a lon & lat map of the upward flux being ejected, overplot the Sub-Earth and Sub-Solar points. 
      release_map = fltarr(360, 180) 
      
      if keyword_set(ICME) then begin
        ; Determine the flux per sq cm surface area if the upward flux were globally uniform .
          average_ejection_globally_integrated = Upward_flux_at_exobase_via_Szabo / (4.*!pi*(re*1.e5)^2)
        
        ; the surface area per bin differs, account for this  
          normalize_area                       = surface_area_per_lat_lon_bin / mean(surface_area_per_lat_lon_bin)
  
          for i = 0, N_Particles - 1 do begin
            release_map[(longitudes[i]*!radeg + 360.) mod 360., lat[i]*!radeg+90.] = release_map[(longitudes[i]*!radeg + 360.) mod 360., lat[i]*!radeg+90.] $
                                                                                     + loc[4,i]*average_ejection_globally_integrated * $
                                                                                     normalize_area[(longitudes[i]*!radeg + 360.) mod 360., lat[i]*!radeg+90.]
          endfor
  
        PRINT, 'Calculated Global Rate in atoms/s = ', TOTAL(surface_area_per_lat_lon_bin*release_map) ; This *MUST* be the same as the global precipitaion x the yield
      
      release_map = rebin(release_map, 60, 30)

      cgPS_Open, filename = directory+'\Surface_reservoir\Release_map.eps', /ENCAPSULATED, xsize = 7.5, ysize = 7.5
      !P.font = 1
      device, SET_FONT = 'Helvetica Bold', /TT_FONT
      !p.charsize = 1.5
  
        axis_format = {xtickinterval:30, ytickinterval:30, ytitle:'Latitude [degrees]', xtitle:'Planetographic W Longitude'}
        cgimage, release_map, /keep_aspect, /axes, xr = [0,360], yr = [-90,90], axkey=axis_format, POSITION=[0.15, 0.50, 0.75, 0.95], minval = min(release_map), maxval = max(release_map)
        cgCOLORBAR, /right, /vertical, title = 'Ejection rate [cm!U-2!N s!U-1!N]', range = minmax(release_map), /fit
        
        cgplot, (subsolar_lon_ET*!radeg + 360.) mod 360., subsolar_lat_ET*!radeg, psym = 16, color = 'red', /overplot         ; the sub-solar point, at simulation end time (not release times)
        cgplot, (subobserver_lon_ET*!radeg + 360.) mod 360., subobserver_lat_ET*!radeg, psym = 16, color = 'blue', /overplot  ; the sub-observer point, at simulation end time (not release times)
        cgtext, (subobserver_lon_ET*!radeg + 360.) mod 360., subobserver_lat_ET*!radeg+10., 'Sub-'+viewpoint, color = 'blue', align = 0.5
        cgtext, (subsolar_lon_ET*!radeg + 360.) mod 360., subsolar_lat_ET*!radeg-16., 'Sub-Solar', color = 'red', align = 0.5
      cgPS_Close

    endif
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
  if ((STRPOS(speed_distribution, 'MBF') eq 0) or (STRPOS(speed_distribution, 'Johnson_2002') eq 0)) then begin
    window, 4
  
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

;---------------------------------CONE-SHAPED EJECTA ANGLES (KNOWN IN METEOR IMPACT DUST)--------------------------------------------------------------------------------- 
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

  ; where in a spherical coordinate system where z is defined as polar_angle_distribution (theta) = 0 is the z axis, phi is the angle relative to the x axis.
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
  if keyword_set(debug) then cgHISTOPLOT, angles_from_normal, binsize = 2., xtitle = 'Ejection Angle from Surface Normal (Degrees)', Ytitle = 'Number of particles'

endif  
return
end