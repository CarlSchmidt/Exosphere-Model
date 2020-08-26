pro radaccel, line, vel, RMAU, wavelength, intensity, ra


;************************************************************************************************************************
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
;    RA:          [CM/S^2] radiation accceleration on a particle due to solar photon resonant scattering 
;
; MODIFICATION HISTORY:
;      C. Schmidt         BU-CSP 4/14/09 Written for Na-D
;      C. Schmidt 2020    Modified for speed, added K-D and Mg 2853A

;*************************************************************************************************************************

  ; Constants
    c = 2.99792458e10   ; speed of light units cm^2/s
    h = 6.62606896e-27  ; Plank's constant
    
  case line of  
    'H-Lyman-a': begin
      m = 1.00794*1.66053886e-24 ;g/atom
      lamda_rest=1215.67
      end
    'O-I' : begin
      m = 15.9994*1.66053886e-24  ;g/atom
      lamda_rest=1304.
      end
    'Na-D' : begin
        m=22.98976928*1.66053886e-24  ;g/atom
        lamda_rest=5893.  ; between the doublet
      end
    'K-D' : begin
       m=39.0983*1.66053886e-24  ;g/atom
       lamda_rest = 7682. ; between the doublet
      end 
    'Mg-2853' : begin
        m=24.305*1.66053886e-24  ;g/atom
        lamda_rest=2852.965
      end
    'H2O' : begin
      ra = 0. / RMAU^2.
      goto, NO_LINE ;No resonant scattering line
      end  
  endcase
  gvalue, line, vel, RMAU, wavelength, intensity, g
  nu_rest = 2.99792458e8/(lamda_rest*1.e-10) ;rest frequency in hz (sodium d is ~508 THz)
  ra = (h*nu_rest*g) / (m*c)
  NO_LINE:
  ;print, 'Radiation acceleration (cm/s^2) =',ra

end