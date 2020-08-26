Pro Call_Generic_Model

;PURSPOSE: Create frames for subsequent animation

COMMON Output_shared, Plot_range, Output_Size_In_Pixels, Output_Title, Center_in_frame, viewpoint, FOV, N_ticks, Tickstep, Observatory, Above_Ecliptic, Boresight_Pixel

;integration_time_array = ((findgen(60)+1)/60.) * 3./24.        ;Integration times for each run (DAYS) 
;Observatory_Array = ['McDonald', 'CASLEO']
Observatory_Array = ['McDonald']
;UTC_Array         = ['Apr 09, 2005 06:19:18', 'Apr 09, 2005 09:19:22']               ; Universal coordinate time to be modeled, time when image taken.
UTC_Array         = ['Apr 09, 2005 06:19:18', 'Apr 09, 2005 09:19:22']   

;How long did all these model runs take to execute?
start_time = systime(/seconds) ;record the time at which the program starts
for j = 0, n_elements(UTC_Array)-1 do begin
  for i = 0, n_elements(Observatory_Array)-1 do begin
    ;Time_range = [0, integration_time_array[i]]
    Observatory_this_run = Observatory_Array[i]
    UTC_this_run = UTC_Array[j]
    Output_title_this_run = Observatory_Array[i]+strmid(UTC_Array[j], 13, 2)
    Generic_Model, Time_range = Time_range, Observatory_this_run = Observatory_this_run, Output_title_this_run = Output_title_this_run, UTC_this_run = UTC_this_run
    print,'**************************************************************************************************************'
    ;print, 'Time Since Simulation Start =',integration_time_array[i] ,'   hours'
  endfor
endfor  
print,'Execution Time =',(systime(/seconds)-start_time)/3600.,'   hours'


end