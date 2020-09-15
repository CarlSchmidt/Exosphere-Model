pro UVVS_Pointing

  ; Restore Tim Cassidy's save file and check it against the PDS MESSENGER UVVS DDRs
  ; Verify that the SPICE geometry calculations agree. 
  ; cf. https://pds-geosciences.wustl.edu/messenger/mess-e_v_h-mascs-3-virs-cdr-caldata-v1/messmas_2001/document/uvvs_cdr_ddr_sis.pdf

  ; Remove existing kernel Data
    CSPICE_KTOTAL, 'all', count
    PRINT, STRCOMPRESS('Cleaning ' + STRING(count) + ' old kernels out of memory . . .')
    i = 0
    WHILE i LT count DO BEGIN
      CSPICE_KDATA, 0, 'all', file, type, source, handle, found
      CSPICE_UNLOAD, file
      i = i + 1
    ENDWHILE
  
  Kernel_Directory = 'C:\Spice\
  CSPICE_FURNSH, STRCOMPRESS(kernel_directory + 'generic_kernels\spk\planets\de431.bsp')              ; SPK (ephemeris kernel) for planets
  CSPICE_FURNSH, STRCOMPRESS(kernel_directory + 'generic_kernels\lsk\naif0010.tls')   
  CSPICE_FURNSH, STRCOMPRESS(kernel_directory + 'MESSENGER\msgr_040803_150430_150430_od431sc_0.bsp')  ; 630 MB !!! SPK Kernal for the mission duration
  CSPICE_FURNSH, STRCOMPRESS(kernel_directory + 'MESSENGER\messenger_2548.tsc')                       ; Most recent spacecraft clock
  CSPICE_FURNSH, STRCOMPRESS(kernel_directory + 'MESSENGER\msgr_v231.tf')                             ; Spacecraft frames & Instrument FOV kernal
  CSPICE_FURNSH, STRCOMPRESS(kernel_directory + 'MESSENGER\msgr_mascs_v100.ti')                       ; Dynamic frames kernal
  CSPICE_FURNSH, STRCOMPRESS(kernel_directory + 'MESSENGER\msgr_1108_v02.bc')                         ; MESSENGER Spacecraft Orientation CK Files, these are monthly: covers DECEMBER 2011 only
  CSPICE_FURNSH, STRCOMPRESS(kernel_directory + 'MESSENGER\msgr_dyn_v600.tf')                         ; MESSENGER dynamic Frames Kernel defining a series of dynamic frames that support data reduction and analysis
  CSPICE_FURNSH, STRCOMPRESS(kernel_directory + 'MESSENGER\pck00010_msgr_v23.tpc')
  
  UTC    = '2011-Aug-04 02:08:37.43'  ;Tim's first datapoint in UTC
  SC     =   -236
  INST   =   -236000
  UVVS   =   -236600
  Body   =   'Mercury'
 
  cspice_bodvrd, Body, 'RADII', 3, Body_radius
  Body_radius = Body_radius[0]

  ; Load UVVS DDR
    UVVS_DDR = read_mascs_ddr('C:\IDL\Generic Model V2\read_write\MESSENGER_UVVS\ud_02_ns_na.dat')   
    UVVS_UTC_TIME = string(UVVS_DDR.UTC_TIME)  ; convert time from byte to string YYDOYTHH:MM:SS.00
  
  ; Load Tim's save file
    restore, 'C:\IDL\Generic Model V2\read_write\MESSENGER_UVVS\orbit 278 radiance for Carl.sav', /verbose
    Tim_utc_time = utc_time
    Tim_RADIANCE_KR = RADIANCE_KR
    
  ; Check for a match between Tim's save file and the PDS  
    start = where(UVVS_UTC_TIME eq Tim_utc_time[0])
    stop = where(UVVS_UTC_TIME eq Tim_utc_time[-1])

    date = Tim_utc_time
    Tim_ET = dblarr(N_elements(Tim_utc_time))
    year = '20' + StrMid(StrTrim(date,2), 0, 2)
    dayofyear = StrMid(StrTrim(date,2), 2, 3)
    time = StrMid(StrTrim(date,2), 6, 11)
    CALDAT, JULDAY(1, dayofyear, year), month, day
    UTC_string = year + '-' + string(month ,format = "(I2.2)") + '-' + string(day ,format = "(I2.2)") + ' ' + time
    for i = 0, N_elements(Tim_utc_time)-1 do begin     
      cspice_utc2et, UTC_string[i], ET
      Tim_ET[i] = ET 
    endfor
    
    date = UVVS_UTC_TIME[Start:Stop]
    PDS_ET = dblarr(N_elements(date))
    year = '20' + StrMid(StrTrim(date,2), 0, 2)
    dayofyear = StrMid(StrTrim(date,2), 2, 3)
    time = StrMid(StrTrim(date,2), 6, 11)
    CALDAT, JULDAY(1, dayofyear, year), month, day
    UTC_string = year + '-' + string(month ,format = "(I2.2)") + '-' + string(day ,format = "(I2.2)") + ' ' + time
    for i = 0, N_elements(date)-1 do begin
      cspice_utc2et, UTC_string[i], ET
      PDS_ET[i] = ET
    endfor
    
    cgplot, Tim_ET, Tim_RADIANCE_KR, color = 'red', psym=4, /xstyle
    cgplot, PDS_ET, UVVS_DDR[Start:Stop].TOTAL_RADIANCE_KR, color = 'blue', psym=5, /overplot
    ;----> looks like a good match between Tim's numbers and what's on the PDS, now let's see if I can calculate the geometry in the PDS from SPICE
        
  ; First verify that Mercury-SC vector can be replicated with SPICE in IAU Mercury coordinates
    print, 'Mercury to MESSENGER vector from the UVVS DDR on the PDS = ', UVVS_DDR[Start].PLANET_SC_VECTOR_TG ; <---replicate this
    cspice_spkpos, 'MESSENGER', PDS_ET, 'IAU_Mercury', 'None', body, ptarg, ltime
    print, 'Mercury to MESSENGER vector via SPICE (IAU Mercury Frame) = ', ptarg[*,0] 
    ; ---> Off by half a km, good enough for government work!

  ; Next, verify the UVVS boresight unit vector
    print, 'Boresight unit vector from the UVVS DDR on the PDS = ', UVVS_DDR[Start].BORESIGHT_UNIT_VECTOR_CENTER_TG  ; <---replicate this
    
    ; get the FOV and frame definition for the instrument
      cspice_getfov, UVVS, 4, shape, frame, bsight, bounds ;Returns bounds in radians UVVS Atmosphere  FOV =  1.0 x 0.04 degrees
    
    ; get the spacecraft clock at this ephemeris time
      SC_Clock = strarr(N_elements(date))
      boresight_vector_IAU_Mercury = dblarr(3, N_elements(date))
      for i = 0, N_elements(date)-1 do begin
        cspice_sce2s, -236, PDS_ET[i], sclkch
        SC_Clock[i] = sclkch
        
        ; cspice_ckgp will require encoded spacecraft clock time, get that encoded time here
        cspice_scencd, -236, SC_Clock[i], sclkdp

        ; Retrieve the 'IAU_Mercury' reference frame to 'INST' reference frame transformation matrix at time sclkdp with a tolerance of 1.e3 sc clock ticks.
        cspice_ckgp, -236000, sclkdp, 1.e3, 'IAU_Mercury', cmat, clkout, found
        cspice_mtxv, cmat, bsight, bore_ref
        boresight_vector_IAU_Mercury[*,i] = bore_ref
      endfor  

      print, 'UVVS Boresight Unit vector (IAU Mercury Frame) = ', boresight_vector_IAU_Mercury[*,0] 
      ;---> A decent match. Crashes if I run it -236600, but it looks like the -236000 spacecraft bus frame & the -236600 UVVS instrument frame are aligned, anyhow, so this works!

  ; Now find the coordinates of the tangent point, verify the altitude of the tangent point and its lat, lon
    print, 'PDS UVVS DDR Boresight tangent altitude [center, min alt, max alt]:', UVVS_DDR[Start].TARGET_ALTITUDE   ; <---replicate this [center, min alt, max alt]
    print, 'PDS UVVS DDR Boresight tangent Planetocentric Longitude:', UVVS_DDR[Start].TARGET_longitude
    print, 'PDS UVVS DDR Boresight tangent Planetocentric Latitude:', UVVS_DDR[Start].TARGET_latitude

    theta = dblarr(N_elements(date))
    SC_to_tangent_point = dblarr(3, N_elements(date))
    for i = 0, N_elements(date)-1 do begin
      SC_to_planet_unit_vector = -ptarg[*,i]/norm(-ptarg[*,i])
      theta[i] = cspice_vsep( boresight_vector_IAU_Mercury[*,i], SC_to_planet_unit_vector) ; angle between the UVVS boresight unit vector and the MESSENGER-planet unit vectors 
      SC_to_tangent_point[*,i] = boresight_vector_IAU_Mercury[*,i] * norm(-ptarg[*,i]) * cos(theta[i])
    endfor
    
    tangent_point_MERC_IAU = ptarg + SC_to_tangent_point

    cspice_reclat, tangent_point_MERC_IAU, radius, lon, lat
    print, 'SPICE Tangent Altitude:', radius[0] - Body_radius
    print, 'SPICE Tangent Planetary Longitude:', lon[0]*!radeg
    print, 'SPICE Tangent Planetary Latitude:', lat[0]*!radeg
stop
end