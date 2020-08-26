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
  CSPICE_FURNSH, STRCOMPRESS(kernel_directory+'generic_kernels\lsk\naif0010.tls')         ; leap seconds kernel
  CSPICE_FURNSH, STRCOMPRESS(kernel_directory+'generic_kernels\pck\pck00010.tpc')         ; Planet rotational states
  CSPICE_FURNSH, STRCOMPRESS(kernel_directory+'generic_kernels\pck\gravity.tpc')          ; GM Gravitational constants
  CSPICE_FURNSH, STRCOMPRESS(kernel_directory+'generic_kernels\spk\planets\de431.bsp')    ; SPK (ephemeris kernel) for planets
  CSPICE_FURNSH, STRCOMPRESS(kernel_directory+'generic_kernels\spk\satellites\sat319.bsp'); SPK (ephemeris kernel) for satellites 
  CSPICE_FURNSH, STRCOMPRESS(kernel_directory+'generic_kernels\fk\heliospheric_v004u.tf') ; Heliospheric *Dynamic* frame, like "GSE" for looking down on solar system 
  CSPICE_FURNSH, STRCOMPRESS(kernel_directory+'Jupiter_System\jup310.bsp')                ; Jupiter isn't in DE431
  
; Load more specialized spice kernels . . . 
  cspice_bodn2c, Body, Body_id, found
  cspice_bods2c, viewpoint, Observer_id, found

; Rosetta Specific Kernels (see aareadme.txt file in kernel directory)
  If Body_id eq 1000012 then begin
    CSPICE_FURNSH, STRCOMPRESS(kernel_directory + 'ROSETTA\kernels\spk\67P_CHURY_GERAS_2004_2016.BSP') ;covers 2003-12-31T23:58:56 to 2015-12-31T23:58:54
    CSPICE_FURNSH, STRCOMPRESS(kernel_directory + 'ROSETTA\kernels\spk\ORHW_______________00122.BSP')  ;Ephemeris data for the Comet Churyumov-Gerasimenko/67P
    CSPICE_FURNSH, STRCOMPRESS(kernel_directory + 'ROSETTA\kernels\spk\ORHR_______________00122.BSP')  ;Rosetta spacecraft predicted and reconstructed cruise ephemeris
    CSPICE_FURNSH, STRCOMPRESS(kernel_directory + 'ROSETTA\kernels\pck\ROS_CGS_RSOC_V03.TPC')      
    CSPICE_FURNSH, STRCOMPRESS(kernel_directory + 'ROSETTA\kernels\pck\ROS_LUTETIA_RSOC_V03.TPC')
    CSPICE_FURNSH, STRCOMPRESS(kernel_directory + 'ROSETTA\kernels\pck\ROS_STEINS_V04.TPC')
    CSPICE_FURNSH, STRCOMPRESS(kernel_directory + 'ROSETTA\kernels\fk\ROS_V19.TF')
  endif

; HST Specific Kernels
  if (Observer_id eq -48) then begin
    CSPICE_FURNSH, STRCOMPRESS(kernel_directory + '1990-01-01_2006-12-31.bsp	')			;SPK (ephemeris kernel) for HST (only through 2006?)
    CSPICE_FURNSH, STRCOMPRESS(kernel_directory + '2006-12-01_2008-05-01.bsp	')			;SPK (ephemeris kernel) for HST (only through 2008?)
  endif

; MESSENGER Specific Kernels
  if (Observer_id eq -236) then begin
    CSPICE_FURNSH, STRCOMPRESS(kernel_directory + 'MESSENGER\msgr_040803_150430_150430_od431sc_0.bsp')  ; 630 MB !!! SPK Kernal for the mission duration
    CSPICE_FURNSH, STRCOMPRESS(kernel_directory + 'MESSENGER\messenger_2548.tsc')                       ; Most recent spacecraft clock
    CSPICE_FURNSH, STRCOMPRESS(kernel_directory + 'MESSENGER\msgr_v231.tf')                             ; Spacecraft frames & Instrument FOV kernal
    CSPICE_FURNSH, STRCOMPRESS(kernel_directory + 'MESSENGER\msgr_mascs_v100.ti')                       ; Dynamic frames kernal
    CSPICE_FURNSH, STRCOMPRESS(kernel_directory + 'MESSENGER\msgr_1108_v02.bc')                         ; MESSENGER Spacecraft Orientation CK Files, these are monthly: covers DECEMBER 2011 only
    CSPICE_FURNSH, STRCOMPRESS(kernel_directory + 'MESSENGER\msgr_dyn_v600.tf')                         ; MESSENGER dynamic Frames Kernel defining a series of dynamic frames that support data reduction and analysis
    CSPICE_FURNSH, STRCOMPRESS(kernel_directory + 'MESSENGER\pck00010_msgr_v23.tpc')
     
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




;STEREO Specific Kernels
if ((Observer_id eq -234) or (Observer_id eq -235)) then begin
  cspice_furnsh, 'C:\ssw\stereo\gen\data\spice\epm\ahead\ahead_2009_050_definitive_predict.epm.bsp' ; Stereo A spk file
  cspice_furnsh, 'C:\ssw\stereo\gen\data\spice\epm\behind\behind_2009_049_definitive_predict.epm.bsp' ; Stereo B spk file
endif
  
CSPICE_KTOTAL, 'all', count
PRINT, STRCOMPRESS('Loaded ' + STRING(count) + ' new SPICE kernel files . . .')
END
