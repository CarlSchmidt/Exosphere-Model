pro write_FK_file

; REVISION HISTORY:
;    Written:                     P. Lierle, Jul 2025

COMMON Output_shared, Plot_range, Output_Size_In_Pixels, Output_Title, Center_in_frame, viewpoint, FOV, N_ticks, Tickstep, Observatory, Above_Ecliptic, Boresight_Pixel, Aperture_Corners
COMMON Model_shared, Body, Ephemeris_time, Obs_Body_Ltime, Parent_ID, Seed, Directory, Particle_data, Line_data, Debug

  OPENW, 1, directory + "loc_prime.tf"
  PRINTF, 1, "LOC_PRIME Frame: Frame kernel specification."
  PRINTF, 1, ""
  PRINTF, 1, "    The +Z axis is aligned with the J2000-referenced Observer to Body vector."
  PRINTF, 1, ""
  PRINTF, 1, "    The component of the Body to Sun vector orthogonal to +Z is aligned with"
  PRINTF, 1, "    the +X axis."
  PRINTF, 1, ""
  PRINTF, 1, "    The +Y axis is the cross product of the +Z axis and the +Y axis."
  PRINTF, 1, ""
  PRINTF, 1, "\begindata"
  PRINTF, 1, ""
  PRINTF, 1, "    FRAME_LOC_PRIME = 1500051"    ; arbitrary; 1400000-2000000 for non NAIF
  PRINTF, 1, "    FRAME_1500051_NAME = 'LOC_PRIME'"
  PRINTF, 1, "    FRAME_1500051_CLASS = 5"
  PRINTF, 1, "    FRAME_1500051_CLASS_ID = 1500051"
  PRINTF, 1, "    FRAME_1500051_CENTER = '" + Body + "'"
  PRINTF, 1, "    FRAME_1500051_RELATIVE = 'J2000'"
  PRINTF, 1, "    FRAME_1500051_DEF_STYLE = 'PARAMETERIZED'"
  PRINTF, 1, "    FRAME_1500051_FAMILY = 'TWO-VECTOR'"
  PRINTF, 1, "    FRAME_1500051_PRI_AXIS = 'Z'"
  PRINTF, 1, "    FRAME_1500051_PRI_VECTOR_DEF = 'OBSERVER_TARGET_POSITION'"
  PRINTF, 1, "    FRAME_1500051_PRI_OBSERVER = '" + viewpoint + "'"
  PRINTF, 1, "    FRAME_1500051_PRI_TARGET = '" + Body + "'"
  PRINTF, 1, "    FRAME_1500051_PRI_ABCORR = 'NONE'"
  PRINTF, 1, "    FRAME_1500051_PRI_FRAME = 'J2000'"
  PRINTF, 1, "    FRAME_1500051_SEC_AXIS = 'X'"
  PRINTF, 1, "    FRAME_1500051_SEC_VECTOR_DEF = 'OBSERVER_TARGET_POSITION'"
  PRINTF, 1, "    FRAME_1500051_SEC_OBSERVER = '" + Body + "'"
  PRINTF, 1, "    FRAME_1500051_SEC_TARGET = 'SUN'"
  PRINTF, 1, "    FRAME_1500051_SEC_ABCORR = 'NONE'"
  PRINTF, 1, "    FRAME_1500051_SEC_FRAME = 'J2000'"
  PRINTF, 1, ""
  WRITEU, 1, "\enddata"   ; WRITEU doesn't append newline char that breaks cspice_furnsh
  CLOSE, 1
  
END