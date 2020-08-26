Pro Call_observing_sequence

  files = [FILE_SEARCH('C:\IDL\Generic Model V2\read_write\Correlation_results\*_Round2.sav'), FILE_SEARCH('C:\IDL\Generic Model V2\read_write\Correlation_results\*_Round3.sav')]
  correlations = fltarr(n_elements(files))
  for l = 0, n_elements(files)-1 do begin
    restore, files[l], /verbose
    correlations[l] = correl
  endfor
  junk = max(correlations, loc)
  print, files[loc]
  print, [strmid(transpose(files[reverse(SORT(correlations))]), 66, 27), transpose(string(correlations[reverse(SORT(correlations))]))]

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
                                
; round 3                                
  Meteor_impact_UTC_array    = ['2011-08-04 01:45:00', '2011-08-04 01:50:00', $
                                '2011-08-04 01:55:00']                        ; index by i
  Plume_temperature_array    = ['3500K']                                      ; index by j
  loop_times_array           = [20]                                           ; index by j
  Surface_distribution_array = ['Point_[300,-60]', 'Point_[330,-60]','Point_[0,-60]', 'Point_[30,-60]', 'Point_[60,-60]', 'Point_[90,-60]', 'Point_[120,-60]', $
                                'Point_[300,-30]', 'Point_[330,-30]','Point_[0,-30]', 'Point_[30,-30]', 'Point_[60,-30]', 'Point_[90,-30]', 'Point_[120,-30]', $
                                'Point_[300,0]', 'Point_[330,0]','Point_[0,0]', 'Point_[30,0]', 'Point_[60,0]', 'Point_[90,0]', 'Point_[120,0]', $
                                'Point_[300,30]', 'Point_[330,30]','Point_[0,30]', 'Point_[30,30]', 'Point_[60,30]', 'Point_[90,30]', 'Point_[120,30]', $
                                'Point_[300,60]', 'Point_[330,60]','Point_[0,60]', 'Point_[30,60]', 'Point_[60,60]', 'Point_[90,60]', 'Point_[120,60]']   ; index by k

  for i = 0, N_elements(Meteor_impact_UTC_array)-1 do begin
    for j = 0, N_elements(Plume_temperature_array)-1 do begin
      for k = 0, N_elements(Surface_distribution_array)-1 do begin

        Observing_Sequence, Meteor_impact_UTC = Meteor_impact_UTC_array[i], Plume_Temperature = Plume_Temperature_array[j], $
                            Surface_distribution = Surface_distribution_array[k], loop_times = loop_times_array[j]
      endfor
    endfor
  endfor

end