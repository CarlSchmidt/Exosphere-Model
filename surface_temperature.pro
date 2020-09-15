Pro Surface_Temperature, Directory, Body, Time, Temperature_map, Silent = Silent

;PURPOSE: Generate a 3D map of the surface temperatures [SS lon, SS lat, Orbital lon] 
;         include thermal intertia, eclipses, internal, heatflow etc.
;
;INPUTS:  Directory = path where files are written.  
;         Body      = the name of the planetary body (string)
;         Time      = UTC time string (This just sets the phase, a datacube of a full diurnal period is output)
;
;OUTPUTS: Temperature map in Kelvin as a datacube and netcdf with dimensions [SS lon, SS lat, Orbital lon] 
;
;CALLING SEQUENCE: e.g., Surface_Temperature, 'Europa', '2016 March 1 00:00', Temperature_map
;
;DEPENDENCIES: SPICE installation, John Spencer's thermprojrs package
;
;WRITTEN: C. Schmidt (LATMOS) 2016 
;
;OPTIONS: Max_Temp keyord also outputs a map the max surface temperature, as a function of geographic longitude, only tested for Mercury.

;---------------------Constants------------------------------------ 
CASE 1 OF
body eq 'Io': begin
    Rotation_rate = 3.5541 ;days SYNODIC DAY
    diurnal_period = 3.5541 ;days SYONDIC DAY
    heatflow = 2.5 * 1.e7 / 1.e4 ;W/m-2 to ergs cm^-2 s-1 from Veeder et al. 1994. 
end ; Io
body eq 'Europa': begin
    ntinc = 800 ;
    Albedo = [ .65, .45 ] ;leading (Jovian Dawn), trailing (Jovian Dusk) approx avg hemispheric albedos in fig 1 of Rathbun et al. 2014
                          ;Note that Spencer et al., 1999 uses 0.55 as bolometric albedo
    Emissivity = .96 ;Emissivity of ice, Marshall et al., 2011
    Thermal_conductivity = 2.52e2 ;ergs cm ^(-1) Kelvin ^(-1) sec ^(-1) ;from Abramov and Spencer, 2008
    Density = 0.92 ;density of Europa ice in g cm^-3 from Abramov and Spencer, 2008
    Specific_Heat = 1.96e7 ;erg g^(-1) K^(-1) ; Vance and Goodman chapter in the Europa book (Univ of Arizona Press)   
    Rotation_rate = 3.5541 ;days SYNODIC DAY
    diurnal_period = 3.5541 ;days SYONDIC DAY
    Thermal_inertia = 7.e4  ;erg cm^-2 s^-1 K^-0.5, Spencer et al. 1999  
    heatflow = 0. ;.05 * 1.e7 / 1.e4 ;W/m-2 upper limit from Spencer et al., 1999 citing G. W. Ojakangas and D. J. Stevenson, Icarus 81, 220 (1989), to ergs cm^-2 s-1 
    ice = 3 ; water ice sublimation cooling
    corrfac = 0.5
    min_dark_temp = 70. ;minimum global temperature, MUST be assumed where thermal balance breaks down near permenant shadow. The lowest nightime temp in Rathbun et al., 2010
    ;'2015-FEB-23 12:01:30'  ;time when Europa is mid Jovian eclipse in its orbit
    ;'2015 FEB 25 06:40:27'  ;time when Europa is at Jovian noon in its orbit. This is orbital longitude zero in Francois' convention
end ; Europa
body eq 'Ganymede': begin
    Albedo = [ .43, .37 ] ; Leading (Jovian Dawn), trailing (Jovian Dusk) Visual albedos in Spencer, Icarus 1987     
                          ; Note that Spencer 1987 thesis uses 0.30 to 0.32 and Buratti (1995) citing Buratti and Ververka (1983) uses a Bond albedo of 0.35
    ;Thermal_inertia = 7e4  ;erg cm^-2 s^-1 K^-0.5, Pappalardo et al. 2004. Jupiter Book. Same as Europa Spencer, PhD thesis and Spencer et al. 1999
    ;Density = 0.92 ;density of Europa ice in g cm^-3 from Abramov and Spencer, 2008 for Europa
    ;ntinc = 650 ; number of time increments in the day. 
    ntinc = 30000 ; number of time increments in the day. 
    Density = [make_array(7, value = .15), make_array(25, value = .92)] ;density in g cm^-3 from Morrison & Chruikshank, 1973 and Abramov and Spencer, 2008 for Europa, respectively
    Thermal_inertia = [make_array(7, value = 2.2e4), make_array(25, value = 7.0e4)] ;Spencer thesis values for Ganymede.
    zarr=[0.1*(1+indgen(7)), .7+0.2*(1+indgen(25))]  ;thickness in cm Spencer thesis gives heat capacity / unit area = 1.6e6 = dens * cp * thickness. 
    Emissivity = .96 ;Emissivity of ice, Marshall et al., 2011
    Thermal_conductivity = 2.52e2 ;ergs cm ^(-1) Kelvin ^(-1) sec ^(-1) ;from Abramov and Spencer, 2008 for Europa
    Specific_Heat = 1.e7 ;erg g^(-1) K^(-1) ; Spencer thesis,  NOTE: 1.96e7 is from the Vance and Goodman chapter in the Europa book (Univ of Arizona Press)
    Rotation_rate = 7.1664  ;days SYNODIC
    diurnal_period = 7.1664 ;days SYNODIC
    ice = 3 ; water ice sublimation cooling
    corrfac = 0.5
    heatflow = 0.;.05 * 1.e7 / 1.e4 ;W/m-2 upper limit from Spencer et al., 1999 citing G. W. Ojakangas and D. J. Stevenson, Icarus 81, 220 (1989), to ergs cm^-2 s-1 
    eclipse = [0., 0.] ;frequency of eclipse in fractional Ganymede days
    min_dark_temp = 80.
    ;'2015-Feb-27 12:24' ;time when Ganymede is mid Jovian eclipse in its orbit
    ;'2015 MAR 03 02:24' ;time when Ganymede is at Jovian noon in its orbit. This is orbital longitude zero in Francois' convention
end ; Ganymede
body eq 'Moon': begin
    ;*****BEWARE: planetary longitude decreases over time, opposite the convention on other bodies!!! 
    Albedo = 0.11 ; Bond albedo, NASA Fact Sheet
    Emissivity = .95 ;Hale & Hapke, 2002.
    ntinc = 10000
    Thermal_conductivity = 200. ;ergs cm ^(-1) Kelvin ^(-1) sec ^(-1) (Glundlach & Blum, 2013; pg. 484) 
    Density = 1.3 ;g cm^(-3) Cremers et al. (1971) and Cremers and Birkebak (1971). Note 1.3 g/cm^-3 also = 3040 kg/m^3 grain density (Warren, 2001) * 0.43 filling factor (Glundlach & Blum 2013)  
    Specific_Heat = 7.11e6 ;erg g^(-1) K^(-1) (Glundlach & Blum 2013), Note: Hemingway and Robie (1973) give 7.6e6 for lunar regolith  
    Rotation_rate = 29.530587981  ;days (synodic)
    diurnal_period = 29.530587981 ;days (synodic)
    Thermal_inertia =  4.3e4 ;erg cm^-2 s^-1 K^-0.5, (Glundlach & Blum 2013) citing Wesselink (1948)
    heatflow = 0.
    min_dark_temp = 70. ;minimum global temperature, MUST be assumed where thermal balance breaks down near permenant shadow.
    ;'2015-Sep-28 03:28' ;time when Moon is fully in Earth umbral eclipse in its orbit
    ;'2016 Sep 01 09:01' ;Solar eclipse / New moon. This is orbital longitude zero in Francois' convention 
end ; Moon
body eq 'Mercury': begin
    Albedo = .106 ;Allen AAQ ;or use Albedo=.142 NASA's Mercury fact sheet
    Emissivity = .95 ;Hale & Hapke, 2002.
    Thermal_conductivity = 1.8e4 ;ergs cm^(-1) Kelvin^(-1) sec^(-1) Wang and Ip, 2008 1.8e4 
    Density = 2.8 ;g cm^(-3) SOURCE??!! Vasavada et al 1999 claim 1.3 g cm^-3 top layer and 1.8 g cm^-3 bottom layer ---->>> probably lower!
    Specific_Heat = 4.6e7 ;erg g^(-1) K^(-1) Wang and Ip, 2008
    Rotation_rate = 58.6462 ;days
    diurnal_period = 175.939 ;days
    thermal_inertia = sqrt(Thermal_conductivity*Density*Specific_Heat) ;using Hemingway's and Vasavada's values
    heatflow = .020 * 10000000./(10.^(4.)) ;W/m-2 from Vasavada et al 1999, to ergs cm^-2 s-1 
    min_dark_temp = 110. ;minimum global temperature, MUST be assumed where thermal balance breaks down near permenant shadow. 
    ntinc = 3000
    corrfac = 5
    ;'Nov 13, 2006 21:33' ;A random perihelion, Francois' code is in TAA and uses this as zero
    ;Notes: thermal inertias in Wang and IP are an order of magnitude higher than the Ksanfomality et al. 2007 review and some other estimates,
    ;       but these values do reproduce a thermal profile matching both Wang and IP, 2008 and Vasavada et al., 1999. They won't match Peplowski et al, 2012.
    ;       Suspect the maximum temperature figure in Peplowski is a misinterpretation of Vasavada et al. 1999 and is NOT correct!    
end ; Mercury
endcase ; constants

  ;**************************generate a map of temperature in latitude and longitude for the ephemeris time**************************************
  start = systime()
    ;Load generic SPICE kernels
     
      ; Clean any lingering kernels out of memory here:
        cspice_ktotal, 'all', count
        Print, 'Deleting ', strtrim(string(count),2), ' old SPICE kernels from memory'
        i=0
        while i lt count do begin
          cspice_kdata, 0, 'all', file, type, source, handle, found
          cspice_unload, file
          i=i+1
        endwhile
      
      ; Load New Kernels
        CSPICE_FURNSH, STRCOMPRESS('C:\SPICE\generic_kernels\lsk\naif0010.tls')         ; leap seconds kernel
        CSPICE_FURNSH, STRCOMPRESS('C:\SPICE\generic_kernels\pck\pck00010.tpc')         ; Planet rotational states
        CSPICE_FURNSH, STRCOMPRESS('C:\SPICE\Jupiter_System\jup309.bsp')   
        CSPICE_FURNSH, STRCOMPRESS('C:\SPICE\generic_kernels\spk\planets\de421.bsp')    ; SPK (ephemeris kernel) for planets
        CSPICE_FURNSH, STRCOMPRESS('C:\SPICE\generic_kernels\spk\satellites\sat319.bsp'); SPK (ephemeris kernel) for satellites 
        cspice_ktotal, 'all', count
        Print, 'Loaded ', strtrim(string(count),2), ' new Spice kernels'

    Print, 'Generating Surface Thermal Map. . .'   

    ;Set up the grid for temperature (Thermal balance tips the scale at SZA's  = 90 degrees in total darkness):
     ;sunlon = 0
     ;sunlat = 0
  
     ;sunlat = (findgen(10)*18 - 81) / !radeg  ;low resolution
     ;sunlon = (findgen(30)*12 - 180) / !radeg ;low resolution, -180 is midnight on the body
  
     ;sunlat = (findgen(30)*6-87) / !radeg  ;medium resolution
     ;sunlon = (findgen(60)*6-180) / !radeg ;medium resolution
     
     ;sunlon = (findgen(120)*3-180) / !radeg ;high resolution
     ;sunlon = (findgen(180)*2-180) / !radeg ;high resolution
     ;sunlat = (findgen(178)-89) / !radeg   ;no interpolation 
     
     sunlat = (findgen(45)*4-88) / !radeg   ;high resolution
     sunlon = (findgen(360)-180) / !radeg   ;no interpolation ---> NEEDED for accuracy!
     Temperature_map = fltarr(n_elements(sunlon), n_elements(sunlat), ntinc)         
        
    ;define an array of times over the last day on that body, times in the body's rest frame
      cspice_str2et, time, ephemeris_time; convert UTC time to ephemeris time (expressed as the number of ephemeris seconds past J2000)   
      cspice_et2utc, ephemeris_time, "C", 0, utc_current
      time_array_insol = ephemeris_time - reverse((findgen(ntinc) / ntinc) * (24.*3600.*diurnal_period))
      insol = time_array_insol ;Initialize array where the insolation vs time will be saved starting with midnight for each time in the for loop below

    ;search and identify for eclipses during this time window 
      eclipse_lightswitch = make_array(ntinc, /integer, value = 1)
      step    = 60.D ;eclipses lasting > 3 min
      MAXWIN  =  1000 
      TIMFMT  = 'YYYY-MON-DD HR:MN:SC.# ::TDB ::RND'
      TIMLEN  =  41
      cnfine = cspice_celld( 2 )
      cspice_wninsd, time_array_insol[0], time_array_insol[-1], cnfine
      result = cspice_celld( MAXWIN*2)
      Case 1 of 
        ((body eq 'Io') or (body eq 'Europa') or (body eq 'Ganymede') or (body eq 'Callisto')): begin
          cspice_gfoclt, 'FULL', 'Jupiter', 'ellipsoid', 'IAU_Jupiter', body, 'ellipsoid', 'IAU_'+strcompress(body), 'None', 'Sun', step, cnfine, result 
          N_eclipses = cspice_wncard( result ) ;number of eclipses during the diurnal period prior to the date specified
        end
        (body eq 'Moon'): begin
          cspice_gfoclt, 'FULL', 'Earth', 'ellipsoid', 'IAU_Earth', body, 'ellipsoid', 'IAU_'+strcompress(body), 'None', 'Sun', step, cnfine, result 
          N_eclipses = cspice_wncard( result ) ;number of eclipses during the diurnal period prior to the date specified
        end
        (body eq 'Mercury'): begin
          N_eclipses = 0  
        end  
      endcase
      
      cspice_timout, [min(time_array_insol), max(time_array_insol)], TIMFMT, TIMLEN, windowstr 
      print, 'Eclipses searched from : ',   windowstr[0], ' to ',   windowstr[1]
      for i= 0L, (N_eclipses - 1L ) do begin ;for each eclipse (if any) adjust times and print findings
          cspice_wnfetd, result, i, left, right
          
          ;the switch to turn off without sunlight
          eclipse_lightswitch[where(((time_array_insol ge left) and (time_array_insol le right)), /NULL, count)] = 0

          cspice_timout, [left, right], TIMFMT, TIMLEN, timstr 
          print, 'Full shadow eclipse found from (body times): ', timstr[0], ' to ',   timstr[1]
          
          ;Correct these to Earth times
            cspice_spkpos, body, left, 'J2000', 'LT', 'Earth', state, Earth_LT
            left = left + Earth_LT
            cspice_spkpos, body, right, 'J2000', 'LT', 'Earth', state, Earth_LT
            right = right + Earth_LT
            cspice_timout, [left, right], TIMFMT, TIMLEN, timstr 
            print, 'Eclipse visible from Earth: ',   timstr[0], ' to ',   timstr[1], ' (Duration:', float(count)*diurnal_period*24./float(ntinc), ' hours)'
      endfor

      cspice_bodvrd, body, 'RADII', 3, radii
      flat = ( radii[0] - radii[2] ) / radii[0] ;flattening coefficiant if the body isn't sphereical
      cspice_bodvrd, body, 'PM', 3, PM
      Sidereal_period = 360./PM[1] ;flattening coefficiant if the body isn't sphereical

      ;Francois' code works in true anomaly angle for Mercury, not orbital longitude from noon like the others
      if body eq 'Mercury' then begin
        TAA = fltarr(n_elements(time_array_insol))
        Sub_Solar_Lon = fltarr(n_elements(time_array_insol))
        cspice_spkezr, body, ephemeris_time, 'J2000', 'None', 'Sun', peri_state, ltime
        for j = 0, n_elements(TAA)-1 do begin ;for each time  
          cspice_spkezr, body, time_array_insol[j], 'J2000', 'None', 'Sun', state, ltime
          TAA[j] = cspice_vsep(state[0:2], peri_state[0:2]) * !radeg
          theta  = cspice_vsep(state[0:2], state[3:5])
          Vel = cos(theta) * norm(state[3:5]) ;scalar projection of the relative velocity along the line of sight from the Sun
          if vel lt 0. then TAA[j] = 360. - TAA[j]
          cspice_subslr, 'Intercept: ellipsoid', body, time_array_insol[j], 'IAU_'+strcompress(body), 'None', 'Sun', sub_solar_point_planet_frame, et_minus_lt, body_WRT_Sun
          cspice_recpgr, body, sub_solar_point_planet_frame, radii[0], flat, lon, lat, alt
          Sub_Solar_Lon[j] = lon*!radeg
        endfor
        print, 'Hack: best to use the real definition of planetary true anomaly angle, and arbitrary input times, or better yet a geometry finder for perihelion!' 
        stop
      endif 

    for i_lat = 0, N_elements(sunlat) - 1 do begin ;step through solar latitudes 
        for i_lon = 0, N_elements(sunlon) - 1 do begin ;step through solar longitudes in the body's local time midnight-dawn-noon-dusk-
          
          ;if i_lat lt 22. then continue
          ;if i_lon lt 180. then continue ;Just for Bob, remove this later! CS 2/28/17
           ;Compute the insolation at this local point, that is sunlat[i_lat], sunlon[i_lon], over it's last diurnal period   

                ;get the geographic point in the body-fixed frame at the modeled time
                cspice_subslr, 'Intercept: ellipsoid', body, ephemeris_time, 'IAU_'+strcompress(body), 'LT+S', 'Sun', sub_solar_point_planet_frame, et_minus_lt, body_WRT_Sun
                cspice_recpgr, body, sub_solar_point_planet_frame, radii[0], flat, lon, lat, alt
                 
                if body eq 'Moon' then abs_lon = lon - sunlon[i_lon] else abs_lon = lon + sunlon[i_lon] ;geographic longitude of the point we'll model in this loop :HACK
                ;abs_lat = lat + sunlat[i_lat] ;geographic latitude of the point we'll model in this loop DANGER of exceeding 90
                abs_lat = sunlat[i_lat] ;geographic latitude of the point we'll model in this loop DANGER of exceeding 90
                cspice_pgrrec, body, abs_lon, abs_lat, alt, radii[0], flat, local_rectan ;convert this location back into rectangular coodinates in the body-fixed frame
                
                cos_sza_angle = fltarr(ntinc) ;build an array of the time history of cos(solar zenith angles) at this location (sunlat[i_lat], sunlon[i_lon])
                for j = 0, n_elements(time_array_insol)-1 do begin ;for each time
                  
                  cspice_spkpos, body, time_array_insol[j], 'J2000', 'None', 'Sun', state, ltime
                  Heliocentric_distance = norm(state) / 149597871. ;AU, Heliocentric distance perhaps changes diurnally, so get it instantaneously
                  insol[j] = 1.374e6 / (Heliocentric_distance^2.) ; solar constant from Hanel et al., Units are erg cm-2 s-1 
                  
                  ;find the SZA, solar zenith angle of the body fixed midnight vector, and this latitude (ASSUMES Sun is 0 latitude)
                  ;cspice_ilumin, "Ellipsoid", body, time_array_insol[j], 'IAU_'+strcompress(body), 'LT+S', 'Earth', local_rectan, TRGEPC, body_WRT_Sun, PHASE, SOLAR_LON, EMISSN 
                  cspice_ilumin, "Ellipsoid", body, time_array_insol[j], 'IAU_'+strcompress(body), 'None', 'Earth', local_rectan, TRGEPC, body_WRT_Sun, PHASE, SOLAR_LON, EMISSN 
                  
                  ;now that the solar flux has been scaled with distance, multiply it by the cos of the solar zenith angle == cos(solar_lon)*cos(solar_lat)
                  cos_sza_angle[j] = cos(SOLAR_LON);*cos(sunlat[i_lat]) ;consistent with the law of cosines in spherical geometries where solar lat and lon defined as (0,0)  
                endfor     
                
                ;shift solar longitude to the point measured 
                insol = insol * cos_sza_angle
                insol[where(cos_sza_angle lt 0., /NULL)] = 0. ;Nightside
                
                ;flip the eclipse lightswitch
                insol = insol*eclipse_lightswitch
                
                ;if needed, interpolate leading and trailing albedos
                if n_elements(albedo) gt 1 then alb = INTERPOL(albedo, [-1., 1.], sin(sunlon[i_lon])) else alb = albedo

                if not keyword_set(silent) then begin
                  window, 0, xs = 600, ys = 400
                  plot, diurnal_period*(findgen(ntinc) / ntinc), insol, $
                    title = 'Insolation at latitude'+strcompress(fix(sunlat[i_lat]*!radeg))+ ' and local time'+ strcompress(string(sunlon[i_lon]*!radeg * 24./360. +12, format = '(F10.1)')) + ' over the last diurnal period', $
                    ytitle = 'erg cm!U-2!N s!U-1!N',xtitle = 'Earth Days', xrange = [0, Diurnal_period], xstyle = 1., charsize = 1.4
                endif

          thermprojrs, tsurf, tod, rhel = Heliocentric_distance, alb = alb, rot = diurnal_period, emvty = Emissivity, cp = Specific_Heat,  $
            rho = density, heatflow = heatflow, nrun=4, ntinc = ntinc, nday = 4, corrfac=corrfac, ice = ice, /blind, insol = insol, zarr = zarr, ti = Thermal_inertia, silent = silent

          if not silent then begin ;make a pretty plot of the temperature history at this geographic location        
            window, 1
            cgplot, (tod*!RADEG)/15. + strcompress(string(sunlon[i_lon]*!radeg * 24./360. + 12)) - 24, tsurf, ytitle = 'Surface temperature at solar latitude'+strcompress(fix(sunlat[i_lat]*!radeg)), $
              xtitle = 'Number of hours ago in '+ body +' local time', xstyle = 1., /ynozero, charsize = 1.5, $   
              title = 'History over diurnal period ending at UT ' + strcompress(STRMID(utc_current, 0, 17)) + ' in ' + body + '''s reference frame' 
          endif
          if total(finite(tsurf)) lt ntinc then stop ;check for NaNs
          Temperature_map[i_lon, i_lat, *] = tsurf ;build the array of surface temps over 1 day [lat, lon, time]           
       endfor
    endfor
    Temperature_map[where(Temperature_map lt min_dark_temp, /NULL)] = min_dark_temp ;Force to everywhere above some threshold. thermal balance needs sunlight, so this is needed.    
    
    ;So far we've been increasing longitude left to right, but this looks funny, so flip it
    Temperature_map = REVERSE(Temperature_map, 1)
    
    ;interpolate it to integer values in degrees longitude and latitude, get the indices where interpolation grid is desired
    interp_lat = interpol(findgen(n_elements(sunlat)), sunlat, (findgen(180) - 90) / !radeg)
    interp_lon = interpol(findgen(n_elements(sunlon)), sunlon, (findgen(360) - 180) / !radeg)
    if body ne 'Mercury' then interp_orbital_lon = findgen(360) * ntinc / 360. else interp_orbital_lon = findgen(ntinc)
    Temperature_map = INTERPOLATE(Temperature_map, interp_lon, interp_lat, interp_orbital_lon, /grid)

    if body eq 'Mercury' then begin
      if not silent then begin ;make a pretty plot of the temerpature history at this geographic location  
        min_t = min(temperature_map)
        max_t = max(temperature_map)
        window, 0, xs = 360, ys = 180 ;wanna see a cool movie?!
        for i = 0, ntinc-1 do tv, bytscl( temperature_map[*,*,i], min_t, max_t)
      endif ;silent
      ;shift them to align in sub-solar planetographic longitude, sub-solar point is still always at pixel [180, 90, *]
      for i = 0, n_elements(interp_orbital_lon)-1 do Temperature_map[*,*,i] = shift(Temperature_map[*,*,i], Sub_Solar_Lon[i] - Sub_Solar_Lon[0])
      if not silent then begin ;make a pretty plot of the temerpature history at this geographic location  
        min_t = min(temperature_map)
        max_t = max(temperature_map)
        window, 0, xs = 360, ys = 180 ;wanna see a cool movie?!
        for i = 0, ntinc-1 do tv, bytscl( temperature_map[*,*,i], min_t, max_t)
      endif ;silent
    endif else begin      
      ;shift them so that the sub-solar point is always at pixel [180, 90, *]
      for i = 0, n_elements(interp_orbital_lon)-1 do Temperature_map[*,*,i] = shift(Temperature_map[*,*,i], i)
    endelse

  ;---------------------------Display and save---------------------------------- 
  loadct, 13
  set_plot,'ps'
  device, bits_per_pixel=8, filename = Directory+strcompress(body)+'_temp_0.eps', $
      /portrait, /color, FONT_SIZE=5., XSIZE=6., YSIZE=4., /encapsulated, SCALE=2., /inches
  !p.charthick=3.
  !p.charsize=2.      
  !p.font=1
  DEVICE, SET_FONT='Times', /TT_FONT 
    cgimage, rebin(Temperature_Map[*,*,0], 720, 360), /interp, /AXIS, AXKEYWORDS = {XSTYLE:1, Xthick:3., Xtitle:'Longitude', xtickinterval:179., $ 
        xtickname:string((90*indgen(5))-180), xticklen:-.01, ySTYLE:1, ythick:3., ytitle:'Latitude', ytickinterval:179., $
        ytickv:(90*findgen(3)), ytickname:string((90*findgen(3))-90, format = '(I5)'), yticklen:-.01}, $
        title='Surface Temperature (!Uo!NK)', /scale, POSITION = [0.15, 0.15, 0.85, 0.85], $
        minvalue = min(Temperature_Map), maxvalue = ceil(max(Temperature_Map)), /KEEP_ASPECT_RATIO       
    cgColorbar, Range=[min(Temperature_Map), ceil(max(Temperature_Map))], /FIT
  DEVICE, /CLOSE
  device, bits_per_pixel=8, filename = Directory+strcompress(body)+'_temp_90.eps', $
    /portrait, /color, FONT_SIZE=5., XSIZE=6., YSIZE=4., /encapsulated, SCALE=2., /inches
  DEVICE, SET_FONT='Times', /TT_FONT 
    cgimage, rebin(Temperature_Map[*,*,90], 720, 360), /interp, /AXIS, AXKEYWORDS = {XSTYLE:1, Xthick:3., Xtitle:'Longitude', xtickinterval:179., $ 
        xtickname:string((90*indgen(5))-180), xticklen:-.01, ySTYLE:1, ythick:3., ytitle:'Latitude', ytickinterval:179., $
        ytickv:(90*findgen(3)), ytickname:string((90*findgen(3))-90, format = '(I5)'), yticklen:-.01}, $
        title='Surface Temperature (!Uo!NK)', /scale, POSITION = [0.15, 0.15, 0.85, 0.85], $
        minvalue = min(Temperature_Map), maxvalue = ceil(max(Temperature_Map)), /KEEP_ASPECT_RATIO       
    cgColorbar, Range=[min(Temperature_Map), ceil(max(Temperature_Map))], /FIT
  DEVICE, /CLOSE
  device, bits_per_pixel=8, filename = Directory+strcompress(body)+'_temp_180.eps', $
      /portrait, /color, FONT_SIZE=5., XSIZE=6., YSIZE=4., /encapsulated, SCALE=2., /inches
    DEVICE, SET_FONT='Times', /TT_FONT 
      cgimage, rebin(Temperature_Map[*,*,180], 720, 360), /interp, /AXIS, AXKEYWORDS = {XSTYLE:1, Xthick:3., Xtitle:'Longitude', xtickinterval:179., $ 
          xtickname:string((90*indgen(5))-180), xticklen:-.01, ySTYLE:1, ythick:3., ytitle:'Latitude', ytickinterval:179., $
          ytickv:(90*findgen(3)), ytickname:string((90*findgen(3))-90, format = '(I5)'), yticklen:-.01}, $
          title='Surface Temperature (!Uo!NK)', /scale, POSITION = [0.15, 0.15, 0.85, 0.85], $
          minvalue = min(Temperature_Map), maxvalue = ceil(max(Temperature_Map)), /KEEP_ASPECT_RATIO       
      cgColorbar, Range=[min(Temperature_Map), ceil(max(Temperature_Map))], /FIT
  DEVICE, /CLOSE
  device, bits_per_pixel=8, filename = Directory+strcompress(body)+'_temp_270.eps', $
      /portrait, /color, FONT_SIZE=5., XSIZE=6., YSIZE=4., /encapsulated, SCALE=2., /inches
    DEVICE, SET_FONT='Times', /TT_FONT 
      cgimage, rebin(Temperature_Map[*,*,270], 720, 360), /interp, /AXIS, AXKEYWORDS = {XSTYLE:1, Xthick:3., Xtitle:'Longitude', xtickinterval:179., $ 
          xtickname:string((90*indgen(5))-180), xticklen:-.01, ySTYLE:1, ythick:3., ytitle:'Latitude', ytickinterval:179., $
          ytickv:(90*findgen(3)), ytickname:string((90*findgen(3))-90, format = '(I5)'), yticklen:-.01}, $
          title='Surface Temperature (!Uo!NK)', /scale, POSITION = [0.15, 0.15, 0.85, 0.85] , $
          minvalue = min(Temperature_Map), maxvalue = ceil(max(Temperature_Map)), /KEEP_ASPECT_RATIO   
      cgColorbar, Range=[min(Temperature_Map), ceil(max(Temperature_Map))], /FIT
  DEVICE, /CLOSE  
  device, bits_per_pixel=8, filename = Directory + strcompress(body) + '_Equatorial.eps', $
     /portrait, /color, FONT_SIZE=5., XSIZE=6., YSIZE=4., /encapsulated, SCALE=2., /inches
     DEVICE, SET_FONT='Times', /TT_FONT 
     cgplot, findgen(360)/15., temperature_map[180,90,*], ytitle = strcompress(body) + ' Equitorial Surface Temperature', xtitle = strcompress('UT ' + STRMID(utc_current, 0, 17) + ' - ' + body +' Local Time'), $
         xrange = [0,24], xstyle = 1. 
  DEVICE, /CLOSE
  SET_PLOT, 'WIN

;format for netcdf    
  if body ne 'Mercury' then begin
    Temperature_Map = congrid(Temperature_Map, 360, 180, 360)     
    orbital_lon = findgen(360)
    n_orbital_lon = n_elements(orbital_lon)
  endif else begin   
    orbital_lon = fltarr(ntinc)
    n_orbital_lon = n_elements(orbital_lon)
  endelse    
  lat = findgen(180) - 90.
  lon = findgen(360) - 180.
  lon = reverse(lon) ;longitude increase with time as the body rotates, so it decreases from left to right
  n_lat = n_elements(lat)   
  n_lon = n_elements(lon)   
  n_temp = n_elements(Temperature_Map)   
  n_TAA = n_elements(TAA)
  n_Sub_Solar_Lon = n_elements(Sub_Solar_Lon)  
;------------------------------------------------
; open file (clobber means it WILL write over an existing file)
;------------------------------------------------
ncid=ncdf_create(directory+body+'_Temperature.nc', /clobber)
;------------------------------------------------
; Define dimensions
;------------------------------------------------
lonid =             ncdf_dimdef(ncid, 'lon', n_lon)
latid =             ncdf_dimdef(ncid, 'lat', n_lat)
orbital_lon_id =    ncdf_dimdef(ncid, 'orbital_lon', n_orbital_lon)
if body eq 'Mercury' then TAA_id =            ncdf_dimdef(ncid, 'True_anomaly', n_TAA) ;Mercury Only
if body eq 'Mercury' then Sub_Solar_Lon_id =  ncdf_dimdef(ncid, 'Sub_Solar_Lon', n_Sub_Solar_Lon) ;Mercury Only 
;------------------------------------------------
; Define variables
;------------------------------------------------
lon2id =            ncdf_vardef(ncid, 'lon', lonid)
lat2id =            ncdf_vardef(ncid, 'lat', latid)
orbital_lon2id =    ncdf_vardef(ncid, 'orbital_lon', orbital_lon_id)
temp2id =           ncdf_vardef(ncid, 'temp', [lonid, latid, orbital_lon_id])
if body eq 'Mercury' then TAA2id =            ncdf_vardef(ncid, 'true_anomaly', TAA_id);Mercury Only
if body eq 'Mercury' then Sub_Solar_Lon_2id = ncdf_vardef(ncid, 'Sub_Solar_Lon', Sub_Solar_Lon_id);Mercury Only
;------------------------------------------------
; Put in attributes 
;------------------------------------------------
ncdf_attput,ncid,lat2id,'title', 'Sub-Solar latitude'
ncdf_attput,ncid,lat2id,'long_name', 'Sub-Solar latitude'
ncdf_attput,ncid,lat2id,'units', 'degrees'
ncdf_attput,ncid,lon2id,'title', 'Sub-Solar longitude'
ncdf_attput,ncid,lon2id,'long_name', 'Sub-Solar longitude'
ncdf_attput,ncid,lon2id,'units', 'degrees'
ncdf_attput,ncid,orbital_lon2id,'title','Orbital longitude'
ncdf_attput,ncid,orbital_lon2id,'long_name','Orbital longitude'
ncdf_attput,ncid,orbital_lon2id,'units','degrees'
ncdf_attput,ncid,temp2id,'title','Surface Temperature'
ncdf_attput,ncid,temp2id,'long_name','Surface Temperature'
ncdf_attput,ncid,temp2id,'units','Kelvin'
if body eq 'Mercury' then ncdf_attput,ncid,TAA2id,'title','True Anomaly'
if body eq 'Mercury' then ncdf_attput,ncid,TAA2id,'long_name','True Anomaly'
if body eq 'Mercury' then ncdf_attput,ncid,TAA2id,'units','degrees'
if body eq 'Mercury' then ncdf_attput,ncid,Sub_Solar_Lon_2id,'title','Sub-Solar planetary longitude'
if body eq 'Mercury' then ncdf_attput,ncid,Sub_Solar_Lon_2id,'long_name','Sub-Solar planetary longitude'
if body eq 'Mercury' then ncdf_attput,ncid,Sub_Solar_Lon_2id,'units','degrees'
;------------------------------------------------
; End define mode
;------------------------------------------------
ncdf_control, ncid, /endef
;------------------------------------------------
; Write variables
;------------------------------------------------
ncdf_varput, ncid, lon2id, lon
ncdf_varput, ncid, lat2id, lat
ncdf_varput, ncid, orbital_lon2id, orbital_lon
ncdf_varput, ncid, temp2id, Temperature_Map
if body eq 'Mercury' then ncdf_varput, ncid, TAA2id, TAA
if body eq 'Mercury' then ncdf_varput, ncid, Sub_Solar_Lon_2id, Sub_Solar_Lon
;------------------------------------------------
; close netcdf file 
;------------------------------------------------
ncdf_close, ncid

;test read
  ID = NCDF_OPEN( Directory+body+'_Temperature.nc', /NOWRITE)
  tempID = ncdf_varid(ID, 'temp')
  latID = ncdf_varid(ID, 'lat')
  lonID = ncdf_varid(ID, 'lon')
  orbital_lonID = ncdf_varid(ID, 'orbital_lon')
  ;TAAID = ncdf_varid(ID, 'true_anomaly')
  ncdf_varget, ID, tempID, Temperature
  ncdf_varget, ID, latID, lat
  ncdf_varget, ID, lonID, lon ;lon increases with time
  ncdf_varget, ID, orbital_lonID, orbital_lon ;orbital Longitude (zero eq Jovian Noon)
  ;ncdf_varget, ID, TAAID, TAA
  ncdf_close, ID
  
;lat and lon of 0,0 is the subsolar point (center of the Europa_Temperature array), longitude increases with time.
finish = systime()
print, start
print, finish

plot, Temperature[*, 90, 0], /ynozero 
oplot, [fltarr(90), reform(Temperature[180, *, 0]), fltarr(90)]
oplot, 80 + (144-80)*(cos((findgen(360)-180.) / !radeg)^.75>0) ;francois' model

stop
End