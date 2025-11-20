PRO LOAD_SPICE, Kernel_Directory

COMMON Model_shared, Body, Ephemeris_time, Seed, Directory, Particle_data, Line_data, Debug
COMMON Output_shared, Plot_range, Output_Size_In_Pixels, Output_Title, Center_in_frame, viewpoint, FOV, N_ticks, Tickstep, Observatory, Above_Ecliptic

  ; Remove existing kernel Data
    CSPICE_KTOTAL, 'all', count
    PRINT, STRCOMPRESS('Cleaning ' + STRING(count) + ' old kernels out of memory . . .')
    i = 0
    WHILE i LT count DO BEGIN
      CSPICE_KDATA, 0, 'all', file, type, source, handle, found
      CSPICE_UNLOAD, file
      i = i + 1
    ENDWHILE
    
; Load generic kernels that are used universally
  Generic_Kernels = kernel_directory + [$
                    'generic_kernels\lsk\naif0010.tls', $           ; leap seconds kernel
                    'generic_kernels\pck\pck00010.tpc', $           ; Planet rotational states
                    'generic_kernels\pck\gravity.tpc', $            ; GM Gravitational constants
                    'generic_kernels\spk\planets\de431.bsp', $      ; SPK (ephemeris kernel) for planets
                    'generic_kernels\spk\satellites\sat319.bsp', $  ; SPK (ephemeris kernel) for satellites 
                    'generic_kernels\fk\heliospheric_v004u.tf', $   ; Heliospheric *Dynamic* frame, like "GSE" for looking down on solar system 
                    'Jupiter_System\jup310.bsp']                    ; For all things Jupiter System that aren't in DE431
  if !VERSION.OS_FAMILY eq ('unix' or 'MacOS') then Generic_Kernels = repstr(Generic_Kernels, '\', '/')                  
  for i = 0, N_elements(Generic_Kernels)-1 do CSPICE_FURNSH, Generic_Kernels[i]

; Now load more the specialized body-specific and observer-specific kernels . . . 
  if strmid(body, 0, 3) eq '100' then Body_id = long(body) else $   ; For comets, the main level input should be the Body_ID number
      cspice_bodn2c, Body, Body_id, found 
  cspice_bods2c, viewpoint, Observer_id, found

; Comet and Rosetta Specific Kernels (see aareadme.txt file in kernel directory)
  if strmid(body, 0, 3) eq '100' then begin
    Comet_SPKs      = file_search('C:\SPICE\Small_Bodies\*.bsp') ; Hack, need to make this line more unix friendly...
    Rosetta_Kernels = kernel_directory + [$
                      'ROSETTA\kernels\spk\67P_CHURY_GERAS_2004_2016.BSP', $ ; Covers 2003-12-31T23:58:56 to 2015-12-31T23:58:54
                      'ROSETTA\kernels\spk\ORHW_______________00122.BSP', $  ; Ephemeris data for the Comet Churyumov-Gerasimenko/67P
                      'ROSETTA\kernels\spk\ORHR_______________00122.BSP', $  ; Rosetta spacecraft predicted and reconstructed cruise ephemeris
                      'ROSETTA\kernels\pck\ROS_CGS_RSOC_V03.TPC', $
                      'ROSETTA\kernels\pck\ROS_LUTETIA_RSOC_V03.TPC', $
                      'ROSETTA\kernels\pck\ROS_STEINS_V04.TPC', $
                      'ROSETTA\kernels\fk\ROS_V19.TF']
                      
                      
    if !VERSION.OS_FAMILY eq ('unix' or 'MacOS') then Rosetta_Kernels = repstr(Rosetta_Kernels, '\', '/')
    for i = 0, N_elements(Rosetta_Kernels)-1 do CSPICE_FURNSH, Rosetta_Kernels[i]
    for i = 0, N_elements(Comet_SPKs)-1 do CSPICE_FURNSH, Comet_SPKs[i]
  endif
  
; Other Comet-specific Kernels
  if (Observer_id eq -48) then begin
    HST_Kernels = kernel_directory + [$
      '1990-01-01_2006-12-31.bsp', $      ; SPK (ephemeris kernel) for HST (only through 2006?)
      '2006-12-01_2008-05-01.bsp']        ; SPK (ephemeris kernel) for HST (only through 2008?)
    if !VERSION.OS_FAMILY eq ('unix' or 'MacOS') then HST_Kernels = repstr(HST_Kernels, '\', '/')
    for i = 0, N_elements(HST_Kernels)-1 do CSPICE_FURNSH, HST_Kernels[i]
  endif

; HST-specific Kernels
  if (Observer_id eq -48) then begin
      HST_Kernels = kernel_directory + [$
        '1990-01-01_2006-12-31.bsp', $      ; SPK (ephemeris kernel) for HST (only through 2006?)
        '2006-12-01_2008-05-01.bsp']        ; SPK (ephemeris kernel) for HST (only through 2008?)
      if !VERSION.OS_FAMILY eq ('unix' or 'MacOS') then HST_Kernels = repstr(HST_Kernels, '\', '/')
      for i = 0, N_elements(HST_Kernels)-1 do CSPICE_FURNSH, HST_Kernels[i]
  endif

; MESSENGER-specific Kernels
  if (Observer_id eq -236) then begin
    MESSENGER_Kernels = kernel_directory + [$
                        'MESSENGER\msgr_040803_150430_150430_od431sc_0.bsp', $  ; 630 MB !!! SPK Kernel for the mission duration                                                              
                        'MESSENGER\messenger_2548.tsc', $                       ; Most recent spacecraft clock                                                                                
                        'MESSENGER\msgr_v231.tf', $                             ; Spacecraft frames & Instrument FOV kernel                                                                   
                        'MESSENGER\msgr_mascs_v100.ti', $                       ; Dynamic frames kernel                                                                                       
                        'MESSENGER\msgr_1108_v02.bc', $                         ; MESSENGER Spacecraft Orientation CK Files, THESE ARE MONTHLY: COVERS AUGUST 2011 ONLY   
                        'MESSENGER\msgr_1210_v02.bc', $                         ; MESSENGER Spacecraft Orientation CK Files, THESE ARE MONTHLY: COVERS October 2012 ONLY        
                        'MESSENGER\msgr_1304_v02.bc', $                         ; MESSENGER Spacecraft Orientation CK Files, THESE ARE MONTHLY: COVERS April 2013 ONLY        
                        'MESSENGER\msgr_dyn_v600.tf', $                         ; MESSENGER dynamic Frames Kernel defining a series of dynamic frames that support data reduction and analysis
                        'MESSENGER\pck00010_msgr_v23.tpc']                                                                                                                                  
                  
    if !VERSION.OS_FAMILY eq ('unix' or 'MacOS') then MESSENGER_Kernels = repstr(MESSENGER_Kernels, '\', '/')
    for i = 0, N_elements(MESSENGER_Kernels)-1 do CSPICE_FURNSH, MESSENGER_Kernels[i]

  ; the Kernels listed in the page 25 example of https://pds-geosciences.wustl.edu/messenger/mess-e_v_h-mascs-3-virs-cdr-caldata-v1/messmas_2001/document/uvvs_cdr_ddr_sis.pdf are:
  ;    msgr20120224.bc
  ;    msgr20120225.bc
  ;    msgr20120226.bc
  ;    msgr_dyn_v600.tf
  ;    msgr_v210.tf
  ;    msgr_mascs_v100.ti
  ;    naif0010.tls
  ;    pck00009_MSGR_v10.tpc
  ;    messenger_1444.tsc
  ;    msgr_de405_de423s.bsp
  ;    msgr_20040803_20140820_od259sc_0.bsp
  
  ;an example label file in the PDS that I found lists
  ;msgr20120427.bc
  ;msgr20120428.bc
  ;msgr20120429.bc
  ;msgr_dyn_v600.tf
  ;msgr_v231.tf
  ;msgr_mascs_v100.ti
  ;naif0011.tls
  ;pck00010_MSGR_v21.tpc
  ;messenger_2548.tsc
  ;msgr_20040803_20150430_od431sc_2.bsp
  ;mess_usgs_151214.tds
endif

; STEREO-specific Kernels
  if ((Observer_id eq -234) or (Observer_id eq -235)) then begin
    STEREO_Kernels = [$
      'C:\ssw\stereo\gen\data\spice\epm\ahead\ahead_2009_050_definitive_predict.epm.bsp', $ ; Stereo A spk file
      'C:\ssw\stereo\gen\data\spice\epm\behind\behind_2009_049_definitive_predict.epm.bsp'] ; Stereo B spk file
    if !VERSION.OS_FAMILY eq ('unix' or 'MacOS') then STEREO_Kernels = repstr(STEREO_Kernels, '\', '/')
    for i = 0, N_elements(STEREO_Kernels)-1 do CSPICE_FURNSH, STEREO_Kernels[i]
  endif
  
CSPICE_KTOTAL, 'all', count
PRINT, STRCOMPRESS('Loaded ' + STRING(count) + ' new SPICE kernel files . . .')
END
