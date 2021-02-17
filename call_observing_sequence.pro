Pro Call_observing_sequence, Display_Correlation_Results = Display_Correlation_Results

if keyword_set(Display_Correlation_Results) then begin
  ;files = FILE_SEARCH('C:\IDL\Generic Model V2\read_write\Correlation_results\*_Round_4.sav')
  files = FILE_SEARCH('C:\IDL\Generic Model V2\read_write\Correlation_results\*_Event2_Round_0.sav')
  ;files = FILE_SEARCH('C:\IDL\Generic Model V2\read_write\Correlation_results\*Round1.sav')
  Correlations_Na = fltarr(n_elements(files)) & Correlations_Mg = fltarr(n_elements(files))
  for l = 0, n_elements(files)-1 do begin
    restore, files[l], /verbose
    correlations_Na[l] = correl_Na
    ;correlations_Mg[l] = correl_Mg
  endfor
  junk = max(correlations_Na, loc)
  print, [strmid(transpose(files[reverse(SORT(correlations_Na))]), 66, strpos(transpose(files[reverse(SORT(correlations_Na))]),'_corr') - 66), $
          transpose(string(correlations_Na[reverse(SORT(correlations_Na))])), transpose(string(correlations_Mg[reverse(SORT(correlations_Na))]))], FORMAT='(A29,F10.5,F10.5)'
  stop        
endif

; ======================================= Post Bug Fix round 1 correlations ==================================================================
;Meteor_impact_UTC_array    = ['2011-08-04 01:40:00', '2011-08-04 01:50:00']           ; index by i
;Plume_temperature_array    = ['3500K']                                                ; index by j
;loop_times_array           = [20]                                                     ; index by j
;Surface_distribution_array = ['Point_[300,-60]', 'Point_[330,-60]', 'Point_[0,-60]', 'Point_[30,-60]',$  ;lon, lat
;                              'Point_[300,-30]', 'Point_[330,-30]', 'Point_[0,-30]', 'Point_[30,-30]',$
;                              'Point_[300,0]',   'Point_[330,0]',   'Point_[0,0]',   'Point_[30,0]']      ; index by k
;                              
;;augmenting until we're not on an edge in space    
;Meteor_impact_UTC_array    = ['2011-08-04 01:40:00', '2011-08-04 01:50:00', '2011-08-04 02:00:00']      ; index by i
;Plume_temperature_array    = ['2500K', '3500K', '5000K']                                                ; index by j
;loop_times_array           = [50, 25, 12]                                                               ; index by j
;Surface_distribution_array = ['Point_[315,-15]']                                                        ; index by k                            

; =================================== End Post Bug Fix round 1 correlations ==================================================================
; 
; ; ======================================= Post Bug Fix round 2 correlations ==================================================================
;Meteor_impact_UTC_array    = ['2011-08-04 01:50:00', '2011-08-04 02:00:00','2011-08-04 02:10:00','2011-08-04 02:20:00', '2011-08-04 02:30:00']           ; index by i
;Plume_temperature_array    = ['20000K']                                               ; index by j
;loop_times_array           = [5]                                                      ; index by j
;Surface_distribution_array = ['Point_[300,-60]', 'Point_[330,-60]', 'Point_[0,-60]', 'Point_[30,-60]',$  ;lon, lat
;                              'Point_[300,-30]', 'Point_[330,-30]', 'Point_[0,-30]', 'Point_[30,-30]',$
;                              'Point_[300,0]',   'Point_[330,0]',   'Point_[0,0]',   'Point_[30,0]']      ; index by k
;
;Meteor_impact_UTC_array    = ['2011-08-04 02:10:00','2011-08-04 02:20:00']           ; index by i
;Plume_temperature_array    = ['20000K']                                               ; index by j
;loop_times_array           = [10]                                                      ; index by j
;Surface_distribution_array = ['Point_[0,-75]', 'Point_[15,-75]', 'Point_[30,-75]',$  ;lon, lat
;                              'Point_[0,-60]', 'Point_[15,-60]', 'Point_[30,-60]',$ 
;                              'Point_[0,-45]', 'Point_[15,-45]', 'Point_[30,-45]']      ; index by k

; =================================== End Post Bug Fix round 2 correlations ==================================================================

; 
; 
; 
; 
; ========================================== Round 1 top correlations ====================================================================
; 01-40-00---0,-30lat---2000K      0.979055
; 01-50-00---0,-30lat---2000K      0.953119
; 01-40-00---0,-30lat---3500K      0.919412
; 01-40-00---0,-30lat---5000K      0.911880
; 01-50-00---0,-30lat---5000K      0.901882
  
;  Meteor_impact_UTC_array    = ['2011-08-04 01:40:00', '2011-08-04 01:50:00',$
;                                '2011-08-04 02:00:00', '2011-08-04 02:10:00']                  ; index by i
;  Plume_temperature_array    = ['2000K','3500K','5000K']                                       ; index by j
;  loop_times_array           = [40, 20, 10]                                                    ; index by j
;  Surface_distribution_array = ['Point_[0,-30]', 'Point_[30,-30]', 'Point_[60,-30]',$
;                                'Point_[0,0]',   'Point_[30,0]',   'Point_[60,0]',$ 
;                                'Point_[0,30]',  'Point_[30,30]', 'Point_[60,30]']             ; index by k     

;  Meteor_impact_UTC_array    = ['2011-08-04 01:40:00', '2011-08-04 01:50:00']                  ; index by i
;  Plume_temperature_array    = ['3500K']                                                       ; index by j
;  loop_times_array           = [20]                                                            ; index by j
;  Surface_distribution_array = ['Point_[345,-45]', 'Point_[0,-45]', 'Point_[15,-45]',$
;                                'Point_[345,-30]', 'Point_[0,-30]', 'Point_[15,-30]',$
;                                'Point_[345,-15]', 'Point_[0,-15]', 'Point_[15,-15]']          ; index by k

;                              01-40-00---0,-50lat---3500K      0.936900
;                              01-50-00---0,-45lat---3500K      0.815195
;                              01-50-00---0,-50lat---3500K      0.813482
;                              01-40-00---15,-45lat---3500      0.775894
;                              01-50-00---15,-45lat---3500      0.752677
;                              01-40-00---345,-30lat---350      0.695242
;                              01-40-00---0,-30lat---3500K      0.557929
;                              01-50-00---15,-15lat---3500      0.541231
;                              01-50-00---345,-30lat---350      0.517642
;                              01-50-00---15,-30lat---3500      0.495098
;                              01-50-00---0,-30lat---3500K      0.415101
;                              01-40-00---15,-30lat---3500      0.411298
;                              01-40-00---345,-15lat---350      0.374178
;                              01-50-00---345,-15lat---350      0.332587
;                              01-40-00---15,-15lat---3500      0.252216

;   Meteor_impact_UTC_array    = ['2011-08-04 01:45:00', '2011-08-04 01:50:00',$
;                                '2011-08-04 01:55:00']                                  ; index by i
;   Plume_temperature_array    = ['3500K']                                      ; index by j
;   loop_times_array           = [20]                                                ; index by j
;   Surface_distribution_array = ['Point_[345,-60]', 'Point_[0,-60]', 'Point_[15,-60]',$
;                                'Point_[345,-45]', 'Point_[0,-45]', 'Point_[15,-45]',$ 
;                                'Point_[345,-30]', 'Point_[0,-30]', 'Point_[15,-30]']   ; index by k   
                                
;; round 3                                
;  Meteor_impact_UTC_array    = ['2011-08-04 01:45:00', '2011-08-04 01:50:00', $
;                                '2011-08-04 01:55:00']                        ; index by i
;  Plume_temperature_array    = ['3500K']                                      ; index by j
;  loop_times_array           = [20]                                           ; index by j
;  Surface_distribution_array = ['Point_[300,-60]', 'Point_[330,-60]','Point_[0,-60]', 'Point_[30,-60]', 'Point_[60,-60]', 'Point_[90,-60]', 'Point_[120,-60]', $
;                                'Point_[300,-30]', 'Point_[330,-30]','Point_[0,-30]', 'Point_[30,-30]', 'Point_[60,-30]', 'Point_[90,-30]', 'Point_[120,-30]', $
;                                'Point_[300,0]', 'Point_[330,0]','Point_[0,0]', 'Point_[30,0]', 'Point_[60,0]', 'Point_[90,0]', 'Point_[120,0]', $
;                                'Point_[300,30]', 'Point_[330,30]','Point_[0,30]', 'Point_[30,30]', 'Point_[60,30]', 'Point_[90,30]', 'Point_[120,30]', $
;                                'Point_[300,60]', 'Point_[330,60]','Point_[0,60]', 'Point_[30,60]', 'Point_[60,60]', 'Point_[90,60]', 'Point_[120,60]']   ; index by k

;------------------------------long hiatus possible bug fixes--------------------------------------------
;; round 4 (runtime = 50 hours)
;  Meteor_impact_UTC_array    = ['2011-08-04 01:35:00', '2011-08-04 01:45:00', $
;                                '2011-08-04 01:55:00', '2011-08-04 02:05:00'] ; index by i
;  Plume_temperature_array    = ['3500K']                                      ; index by j
;  loop_times_array           = [40]                                           ; index by j
;  Surface_distribution_array = ['Point_[0,0]', 'Point_[30,0]','Point_[60,0]', 'Point_[90,0]', 'Point_[120,0]', 'Point_[150,0]', 'Point_[180,0]', $
;                                'Point_[210,0]', 'Point_[240,0]','Point_[270,0]', 'Point_[300,0]', 'Point_[330,0]']   ; index by k

;; round 5 (runtime = ?? hours)
;  Meteor_impact_UTC_array    = ['2011-08-04 02:10:00', $
;                                '2011-08-04 02:15:00', '2011-08-04 02:20:00'] ; index by i
;  Plume_temperature_array    = ['2000K', '3500K', '5000K']                    ; index by j
;  loop_times_array           = [150, 80, 20]                                  ; index by j
;  Surface_distribution_array = ['Point_[100,-15]', 'Point_[115,-15]', $
;                                'Point_[100,0]', 'Point_[115,0]', $
;                                'Point_[100,15]', 'Point_[115,15]' ]   ; index by k


;; round 6 (runtime = ??? Full resolution)
;  Meteor_impact_UTC_array    = ['2011-08-04 02:15:00']                         ; index by i
;  Plume_temperature_array    = ['3502K', '5000K']                              ; index by j
;  loop_times_array           = [100, 80]                                       ; index by j
;  Surface_distribution_array = ['Point_[100,0]', 'Point_[115,0]', 'Point_[100,15]', 'Point_[115,15]' ]   ; index by k

; Event 2 round 0 (runtime = ??? Full resolution)
;  Meteor_impact_UTC_array    = ['2013-04-13 13:15:00', '2013-04-13 13:25:00', '2013-04-13 13:35:00', '2013-04-13 13:45:00']                         ; index by i
;  Plume_temperature_array    = ['10000K']                              ; index by j
;  loop_times_array           = [70]                                    ; index by j
;  Surface_distribution_array = ['Point_[150, 0]', 'Point_[180, 0]', 'Point_[210, 0]', 'Point_[240, 0]', 'Point_[270, 0]' ]   ; index by k

;  Meteor_impact_UTC_array    = ['2013-04-13 13:35:00', '2013-04-13 13:40:00', '2013-04-13 13:45:00']                         ; index by i
;  Plume_temperature_array    = ['10000K']                              ; index by j
;  loop_times_array           = [20]                                    ; index by j
;  Surface_distribution_array = ['Point_[255, 0]', 'Point_[270, 0]', 'Point_[270, 20]' ]   ; index by k

  Meteor_impact_UTC_array    = ['2013-04-13 13:35:00']                         ; index by i
  Plume_temperature_array    = ['10000K', '15000K']                              ; index by j
  loop_times_array           = [40,40]                                    ; index by j
  Surface_distribution_array = ['Point_[235, 50]','Point_[235, 30]','Point_[245, 30]','Point_[225, 30]']   ; index by k
  for i = 0, N_elements(Meteor_impact_UTC_array)-1 do begin
    for j = 0, N_elements(Plume_temperature_array)-1 do begin
      for k = 0, N_elements(Surface_distribution_array)-1 do begin
        ;if ((i eq 1) and (k lt 4)) then continue
        Observing_Sequence_Event_3, Meteor_impact_UTC = Meteor_impact_UTC_array[i], Plume_Temperature = Plume_Temperature_array[j], $
                                    Surface_distribution = Surface_distribution_array[k], loop_times = loop_times_array[j]
      endfor
    endfor
  endfor
Observing_Sequence_Event_2

end