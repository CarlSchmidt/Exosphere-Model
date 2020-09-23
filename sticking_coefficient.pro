FUNCTION STICKING_COEFFICIENT, species, temperature

  ; The following data was digitized from Yakshinskiy and Madey, (2005), https://doi.org/10.1016/j.susc.2005.06.062, Fig. 3
  ; Sticking probablilities for Na and K from a SiO2 substrate of different temperatures
    Case Species of
      'K':  S_vs_T = [[99.89465859571243, 1.0028008911736639],$
                      [300.5225309499101, 0.9781779382317259],$
                      [399.8086315222407, 0.8417582096346145],$
                      [500.5177284676532, 0.969187691446779]]
      'Na': S_vs_T = [[99.89465859571243, 1.0028008911736639],$
                      [249.62561518753697, 0.49915163106913685],$
                      [300.20097344227315, 0.3762222839352791],$
                      [500.10742944005153, 0.20110791177631282]]
      'Mg': S_vs_T = [[99.89465859571243, 1.0028008911736639],$   ; HACK HACK!!! ASSUME MG BOUNCES LIKE NA, NOT AT ALL A JUSTIFIED ASSUMPTION HACK HACK!!! 
                      [249.62561518753697, 0.49915163106913685],$
                      [300.20097344227315, 0.3762222839352791],$
                      [500.10742944005153, 0.20110791177631282]]                
                      
       Else: S_vs_T = [[10., 1.],$                                ; For other species I'm unsure, easiest to effectively switch off bouncing
                      [1000., 1.],$
                      [10000., 1.]]            
    endcase 
    
  ; Generate a fit to the sticking probablitiy.
    T = S_vs_T[0,*]                                 ; temperature data
    S = S_vs_T[1,*]                                 ; sticking probablility data
    coefficients  = POLY_FIT(T, S, 2.)
    sticking_prob = POLY(temperature, coefficients)
    
  ; If there's no data at the input temperature(s), don't extrapolate
    sticking_prob[where(temperature gt max(T), /NULL)] = min(S)
    sticking_prob[where(temperature lt min(T), /NULL)] = 1. 

RETURN, Sticking_prob
END
