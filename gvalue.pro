pro GVALUE, line, vel, RMAU, wavelength, intensity, g



; Check Killen et al. (2009)

;RMAU = .352
;Vel  = -12000.


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
;
;!!!!O and H yeild values which differ from those published by as much as 1/20 and 2 for each respectively.

  case line of  
    'H-Lyman-a': begin
      spectrum_portion='spectre_1210_10.dat' ;H Lyman Alpha 
      f1=.2776    ;1215.668 A oscillator stretgths
      f2=.1388    ;1215.674 A
      f=f1+f2 ;unitless total oscillator strength for H Lyman Alpha transistion (from NIST website)
      lamda_rest=1215.67
      fudge=3.75e-6 ;fudge factor to match Killen et al., 08!!!!!
      end
    'O-I': begin
      spectrum_portion='spectre_1300_10.dat'  ;OI @ 1304A
      f1=.0520   ;1302.168 A oscillator stretgths
      f2=.0518   ;1304.858 A
      f3=.0519   ;1306.029 A
      f=f1+f2+f3
      lamda_rest=1304.
      fudge=5.3e-4 ;fudge factor to match Killen et al., 08!!!!!
      end
    'Na-D': begin
      lamda_rest2 = 5891.583264 ; Vacuum D2 wavelength from NIST ASD
      lamda_rest1 = 5897.558147 ; Vacuum D1 wavelength from NIST ASD
      f1          = .320        ; Unitless total oscillator strength for D1 5896 transition (from NIST ASD website)
      f2          = .641        ; Unitless total oscillator strength for D2 5890 transition (from NIST ASD website)
      end  
    'K-D' : begin
      lamda_rest2 = 7667.008906 ; Vacuum D2 wavelength from NIST ASD
      lamda_rest1 = 7701.083536 ; Vacuum D1 wavelength from NIST ASD
      f1          = 0.3201      ; Unitless total oscillator strength for D1 7699 transition (from NIST ASD website)
      f2          = 0.6661      ; Unitless total oscillator strength for D2 7665 transition (from NIST ASD website)
      end  
    'Mg-2853' : begin
        lamda_rest = 2852.965   ; Vacuum wavelength from NIST ASD
        f          = 1.80       ; Unitless total oscillator strength for Mg 1 ground state resonance transition (from NIST ASD website)
      end  
  endcase

if (line eq 'H2O') then begin ;No resonant scattering line, returns 1.d scattering --> emission will be 1.e-6 * the cm column density if in Rayleighs
  g = replicate(1.d, size(vel, /dimensions))
endif

; ------------------------------Constants-----------------------------------------
  h      = 6.62606876e-27           ; Plank's constant in ergs*s
  e      = 4.80320427e-10           ; Electron charge in esu or statcoulombs
  me     = 9.10938188e-28           ; Electron mass g 
  c      = 2.99792458e10            ; Speed of light in cm/s
  esigma = (!pi*(e^2.))/(me*c)      ; Classical electron cross-section in cm^2*Hz/atom
  c      = 2.99792458e8             ; speed of light units m^2/s
  h      = 6.62606896e-34 

if (line eq 'Na-D') or (line eq 'K-D') then begin  ;Sodium has adjacent 2 spectral lines, so they can be calculated as a pair  
  nu_rest1              = c / (double(lamda_rest1)*1.e-10) ;rest frequency in vacuum 
  nu_rest2              = c / (double(lamda_rest2)*1.e-10) ;rest frequency
  lamda_particle_frame1 = double(lamda_rest1) + double(lamda_rest1)*(vel/c) ; Doppler shifted wavelengths 
  lamda_particle_frame2 = double(lamda_rest2) + double(lamda_rest2)*(vel/c) ; Doppler shifted wavelengths 
  a1                    = f1 * esigma ; absorption cross-section for this transition in m^2*Hz/atom
  a2                    = f2 * esigma ; absorption cross-section for this transition in m^2*Hz/atom
  g1                    = interpol(intensity, wavelength, lamda_particle_frame1) * a1 * double(lamda_rest1) / ( nu_rest1 * RMAU^2 ) ; g-value in photons/atom/s for D1, lamda/v is a conversion factor to match units
  g2                    = interpol(intensity, wavelength, lamda_particle_frame2) * a2 * double(lamda_rest2) / ( nu_rest2 * RMAU^2 ) ; g-value in photons/atom/s for D2, lamda/v is a conversion factor to match units
  g                     = g1 + g2
  ;print, mean(g1), mean(g2)
endif else begin
  nu_rest               = c / (lamda_rest*1.e-10) ;rest frequency
  lamda_particle_frame  = lamda_rest + lamda_rest*(vel/c);Doppler shifted wavelengths  
  a                     = f*esigma
  g                     = interpol(intensity, wavelength, lamda_particle_frame) * a * double(lamda_rest) / ( nu_rest * RMAU^2 ) ; g-value in photons/atom/s for D1, lamda/v is a conversion factor to match units
endelse

return
end