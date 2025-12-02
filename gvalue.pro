pro GVALUE, line, vel, RMAU, wavelength, intensity, g

;+
; NAME:
;     GVALUE
;
; PURPOSE:
;     Compute scattering rates for sunlight in photons per atom per second. 
;
; CALLING SEQUENCE:
;     GVALUE, line, vel, RMAU, wavelength, intensity, g
;
; INPUT PARAMETERS:
;    LINE:        case sensitive string 'H-LYMAN-A', 'O-I', 'NA-D', 'K-D' are accepted 
;    VEL:         [M/S] heliocentric velocity 
;    RMAU:        [AU] distance from the sun to scattering particle, can be an array 
;    WAVELENGTH:  [ANGSTROMS] ***IN VACCUUM*** wavelength array of a local portion of the solar spectrum with correspnding INTENSITY  
;    INTENSITY:   [PHOTONS / (CM^2 * ANGSTROMS * S)] solar spectral irradiance at 1 AU
; RETURNS:
;    G:           [PHOTONS / (ATOM*S)] the mean photon scattering rate
;
; MODIFICATION HISTORY:
;      C. Schmidt         BU-CSP 4/14/09 Written for Na-D
;      C. Schmidt 2020    Modified for speed, added K-D and Mg 2853A
;      C. Schmidt 2021    Modified for spline interpolation, added Ca 4227A
;      C. Schmidt 2022    Added several Ca+, Al, Mn, and the blue K doublet
;
;!!!!O and H yeild values which differ from those published by as much as 1/20 and 2 for each respectively.

SetDefaultValue, w1, 1.
SetDefaultValue, w2, 1.
SetDefaultValue, w3, 1.

  case line of  
;    'H-Lyman-a': begin
;      N_lines = 1
;      spectrum_portion='spectre_1210_10.dat' ;H Lyman Alpha 
;      f1=.2776    ;1215.668 A oscillator stretgths
;      f2=.1388    ;1215.674 A
;      f=f1+f2 ;unitless total oscillator strength for H Lyman Alpha transistion (from NIST website)
;      lamda_rest=1215.67
;      fudge=3.75e-6 ;fudge factor to match Killen et al., 08!!!!!
;      end
;    'O-I': begin
;      N_lines     = 3
;      spectrum_portion='spectre_1300_10.dat'  ;OI @ 1304A
;      f1=.0520   ;1302.168 A oscillator strengths
;      f2=.0518   ;1304.858 A
;      f3=.0519   ;1306.029 A
;      f=f1+f2+f3
;      lamda_rest=1304.
;      fudge=5.3e-4 ;fudge factor to match Killen et al., 08!!!!!
;      end
    'Na-D': begin
      N_lines     = 2
      lamda_rest2 = 5891.583264 ; Vacuum D2 wavelength from NIST ASD
      lamda_rest1 = 5897.558147 ; Vacuum D1 wavelength from NIST ASD
      f1          = .320        ; Unitless total oscillator strength for D1 5896 transition (from NIST ASD website)
      f2          = .641        ; Unitless total oscillator strength for D2 5890 transition (from NIST ASD website)
      end  
    'K-D' : begin
      N_lines     = 2
      lamda_rest2 = 7667.008906 ; Vacuum D2 wavelength from NIST ASD
      lamda_rest1 = 7701.083536 ; Vacuum D1 wavelength from NIST ASD
      f1          = 0.3201      ; Unitless total oscillator strength for D1 7699 transition (from NIST ASD website)
      f2          = 0.6661      ; Unitless total oscillator strength for D2 7665 transition (from NIST ASD website)
      end  
    'K-Blue' : begin
      N_lines     = 2
      lamda_rest2 = 4045.279  ; Vacuum wavelength from NIST ASD
      lamda_rest1 = 4048.351  ; Vacuum wavelength from NIST ASD
      f1          = 5.67e-03  ; Unitless total oscillator strength (from NIST ASD website)
      f2          = 2.63e-03  ; Unitless total oscillator strength (from NIST ASD website)
      end
    'Mg-2853' : begin
      N_lines    = 1
      lamda_rest = 2852.965   ; Vacuum wavelength from NIST ASD
      f          = 1.80       ; Unitless total oscillator strength for Mg I ground state resonance transition (from NIST ASD website)
      end
    'Li-D'  : begin
      N_lines     = 2
      lamda_rest1 = 6707.76 ; Vacuum D2 wavelength from NIST ASD
      lamda_rest2 = 6707.91 ; Vacuum D1 wavelength from NIST ASD
      f1          = 4.9797e-01      ; Unitless total oscillator strength for D transition (from NIST ASD website)
      f2          = 2.4899e-01      ; Unitless total oscillator strength for D transition (from NIST ASD website)
    end
    'Ca-4227' : begin
      N_lines    = 1
      lamda_rest = 4227.918   ; Vacuum wavelength from NIST ASD
      f          = 1.75       ; Unitless total oscillator strength for Ca I ground state resonance transition (from NIST ASD website)
      end 
    'Al' : begin              ; Vervack et al. (2016) measured line ratio is 2:1, as are his reported g-values
      N_lines     = 2
      lamda_rest1 = 3945.1222 ; Vacuum wavelength from NIST ASD
      lamda_rest2 = 3962.6409 ; Vacuum wavelength from NIST ASD
      f1          = 1.16e-01  ; Unitless total oscillator strength (from NIST ASD website)
      f2          = 1.16e-01  ; Unitless total oscillator strength (from NIST ASD website)
      A1          = 4.99e+07  ; Einstein Aki (from NIST ASD website)
      A2          = 9.85e+07  ; Einstein Aki (from NIST ASD website)
      w1          = A1/(A1+A2); Al I has a multiple J term levels in the ground state, so the branching ratio of the Einstein A matters.
      w2          = A2/(A1+A2); See Chamblain & Hunten (1987) book Eqns 6.2.2 & 6.2.4.
      ; No other resonance lines considered herein have split ground states yet, but watch out for the O 1304 triplet!!!!
      end 
    'Mn' : begin
      N_lines     = 3
      lamda_rest1 = 4031.892   ; Vacuum wavelength from NIST ASD (Ritz)
      lamda_rest2 = 4034.202   ; Vacuum wavelength from NIST ASD (Ritz)
      lamda_rest3 = 4035.623   ; Vacuum wavelength from NIST ASD (Ritz) 
      f1          = 5.5e-02   ; Unitless total oscillator strength (from NIST ASD website)
      f2          = 4.03e-02  ; Unitless total oscillator strength (from NIST ASD website)
      f3          = 2.57e-02  ; Unitless total oscillator strength (from NIST ASD website)
      end
    'Ca_H_and_K' : begin      ; Vervack et al. (2016) measured line ratio is 2:1, as are his reported g-values
      N_lines     = 2
      lamda_rest1 = 3934.777  ; Vacuum wavelength from NIST ASD (Ritz)
      lamda_rest2 = 3969.591  ; Vacuum wavelength from NIST ASD (Ritz)
      f1          = 6.82e-01  ; Unitless total oscillator strength (from NIST ASD website)
      f2          = 3.3e-01   ; Unitless total oscillator strength (from NIST ASD website)
    end
    'Mg-II' : begin      ; Vervack et al. (2016) measured line ratio is 2:1, as are his reported g-values
      N_lines     = 2
      lamda_rest1 = 2798.823  ; Vacuum wavelength from NIST ASD (Ritz)
      lamda_rest2 = 2803.530  ; Vacuum wavelength from NIST ASD (Ritz)
      f1          = 6.08e-01  ; Unitless total oscillator strength (from NIST ASD website)
      f2          = 3.03e-01   ; Unitless total oscillator strength (from NIST ASD website)
    end
    'Fe-I' : begin      
      N_lines     = 3
      lamda_rest1 = 3192.5817 ; these vacuum WLs all go from split ground states into (air) wavelength 5167.4881A, tentatively observed in comet 3I
      lamda_rest2 = 3235.5466  
      lamda_rest3 = 3265.9880      
      f1          = 2.82e-04  ; Unitless total oscillator strength (from NIST ASD website)
      f2          = 1.80e-04  ; Unitless total oscillator strength (from NIST ASD website)
      f3          = 1.99e-04  ; Unitless total oscillator strength (from NIST ASD website)
      A1          = 2.37e+05  ; Einstein Aki (from NIST ASD website)
      A2          = 1.15e+05  ; Einstein Aki (from NIST ASD website)
      A3          = 8.89e+04  ; Einstein Aki (from NIST ASD website)
      w1          = A1/(A1+A2+A3) ; Fe I has a multiple J term levels in the ground state, so the branching ratio of the Einstein A matters.
      w2          = A2/(A1+A2+A3) ; See Chamblain & Hunten (1987) book Eqns 6.2.2 & 6.2.4.
      w3          = A3/(A1+A2+A3)

;      N_lines     = 2
;      lamda_rest1 = 3194.1488 ; these vacuum WLs all go from split ground states into (air) wavelength 5171.5960 , tentatively observed in comet 3I
;      lamda_rest2 = 3237.1561
;      f1          = 6.77e-04  ; Unitless total oscillator strength (from NIST ASD website)
;      f2          = 4.40e-04  ; Unitless total oscillator strength (from NIST ASD website)
;      A1          = 4.43e+05  ; Einstein Aki (from NIST ASD website)
;      A2          = 2.18e+05  ; Einstein Aki (from NIST ASD website)
;      w1          = A1/(A1+A2) ; Fe I has a multiple J term levels in the ground state, so the branching ratio of the Einstein A matters.
;      w2          = A2/(A1+A2) ; See Chamblain & Hunten (1987) book Eqns 6.2.2 & 6.2.4.
    end    

    else : print, 'Requested line has undefined parameters'
  endcase

; ------------------------------Constants-----------------------------------------
  h      = 6.62606876e-27           ; Plank's constant in ergs*s
  e      = 4.80320427e-10           ; Electron charge in esu or statcoulombs
  me     = 9.10938188e-28           ; Electron mass g 
  c      = 2.99792458e10            ; Speed of light in cm/s
  esigma = (!pi*(e^2.))/(me*c)      ; Classical electron cross-section in cm^2*Hz/atom
  c      = 2.99792458e8             ; speed of light units m^2/s
  h      = 6.62606896e-34 

Case N_lines of
  1 : begin
    nu_rest               = c / (lamda_rest*1.e-10)                           ; rest frequency (vacuum)
    lamda_particle_frame  = lamda_rest + lamda_rest*(vel/c)                   ; Doppler shifted wavelengths
    a                     = f*esigma                                          ; absorption cross-section for this transition in m^2*Hz/atom
    g                     = interpol(intensity, wavelength, lamda_particle_frame, /spline) * a * double(lamda_rest) / ( nu_rest * RMAU^2 ) * w1 ; g-value in photons/atom/s for D1, lamda/v is a conversion factor to match units
  end
  2 : begin
    nu_rest1              = c / (double(lamda_rest1)*1.e-10)                  ; rest frequency in vacuum
    nu_rest2              = c / (double(lamda_rest2)*1.e-10)                  ; rest frequency in vacuum
    lamda_particle_frame1 = double(lamda_rest1) + double(lamda_rest1)*(vel/c) ; Doppler shifted wavelengths
    lamda_particle_frame2 = double(lamda_rest2) + double(lamda_rest2)*(vel/c) ; Doppler shifted wavelengths
    a1                    = f1 * esigma                                       ; absorption cross-section for this transition in m^2*Hz/atom
    a2                    = f2 * esigma                                       ; absorption cross-section for this transition in m^2*Hz/atom
    g1                    = interpol(intensity, wavelength, lamda_particle_frame1) * a1 * double(lamda_rest1) / ( nu_rest1 * RMAU^2 ) * w1 ; g-value in photons/atom/s, lamda/v is a conversion factor to match units
    g2                    = interpol(intensity, wavelength, lamda_particle_frame2) * a2 * double(lamda_rest2) / ( nu_rest2 * RMAU^2 ) * w2 ; g-value in photons/atom/s, lamda/v is a conversion factor to match units
    g                     = g1 + g2                                           ; sum the doublet's g-value
  end
  3 : begin
    nu_rest1              = c / (double(lamda_rest1)*1.e-10)                  ; rest frequency in vacuum
    nu_rest2              = c / (double(lamda_rest2)*1.e-10)                  ; rest frequency in vacuum
    nu_rest3              = c / (double(lamda_rest3)*1.e-10)                  ; rest frequency in vacuum
    lamda_particle_frame1 = double(lamda_rest1) + double(lamda_rest1)*(vel/c) ; Doppler shifted wavelengths
    lamda_particle_frame2 = double(lamda_rest2) + double(lamda_rest2)*(vel/c) ; Doppler shifted wavelengths
    lamda_particle_frame3 = double(lamda_rest3) + double(lamda_rest3)*(vel/c) ; Doppler shifted wavelengths
    a1                    = f1 * esigma                                       ; absorption cross-section for this transition in m^2*Hz/atom
    a2                    = f2 * esigma                                       ; absorption cross-section for this transition in m^2*Hz/atom
    a3                    = f3 * esigma                                       ; absorption cross-section for this transition in m^2*Hz/atom
    g1                    = interpol(intensity, wavelength, lamda_particle_frame1) * a1 * double(lamda_rest1) / ( nu_rest1 * RMAU^2 ) * w1 ; g-value in photons/atom/s, lamda/v is a conversion factor to match units
    g2                    = interpol(intensity, wavelength, lamda_particle_frame2) * a2 * double(lamda_rest2) / ( nu_rest2 * RMAU^2 ) * w2 ; g-value in photons/atom/s, lamda/v is a conversion factor to match units
    g3                    = interpol(intensity, wavelength, lamda_particle_frame3) * a3 * double(lamda_rest3) / ( nu_rest3 * RMAU^2 ) * w3 ; g-value in photons/atom/s, lamda/v is a conversion factor to match units
    g                     = g1 + g2 + g3                                      ; sum the triplet's g-value
  end
  else : print, 'Undefined Line'
Endcase

return
end