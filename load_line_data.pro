Pro LOAD_LINE_DATA, LINE
  ;+
  ; NAME:
  ;     LOAD_LINE_DATA
  ;     
  ; PURPOSE:
  ;     Load in a high resolution solar spectrum around the line feature of interest. 
  ;     
  ; CALLING SEQUENCE:
  ;     LOAD_LINE_DATA, LINE
  ;
  ; INPUT PARAMETERS:
  ;     LINE    = String
  ;
  ; RETURNS:
  ;     Line_data = Structure consisting of
  ;                 .name, .line, .unity_optical_depth, .wavelength, .intensity
  ;
  ; MODIFICATION HISTORY:
  ;      C. Schmidt 2010    Written, used Killen et al. 2009 ApJs for solar flux
  ;      C. Schmidt 2020    Modified for Kurucz model spectra

  COMMON Model_shared, Body, Ephemeris_time, Seed, Directory, Particle_data, line_data, Debug
  
  ;===================================FORMAT NAME STRINGS=================================================
  
  CASE Line OF 
    'H-Lyman-a': NAME = strcompress('H Lyman '+cgSymbol('alpha'))
    'OI':        NAME = strcompress('Oxygen 1304'+ cgSymbol('angstrom'))
    'Na-D':      NAME = 'Na D!D1!N + D!D2!N'
    'K-D':       NAME = 'K D!D1!N + D!D2!N'
    'Mg-2853':   NAME = strcompress('MG 2853'+ cgSymbol('angstrom'))
    'H2O':       NAME = 'H!D2!NO' 
  ENDCASE
  
  ;==================================UNITY OPTICAL DEPTH===================================================
  
  ;NOTE: These are rough numbers, See Killen et al., 2009
  CASE Line OF  ;column density for an optical depth of unity cm^-2
    'H-Lyman-a': unity_optical_depth = 3.6077e12
    'OI':        unity_optical_depth = 1.7424e13
    'Na-D':      unity_optical_depth = 1.645e11
    'K-D':       unity_optical_depth = 9.e99 ; TBD
    'Mg-2853':   unity_optical_depth = 9.e99 ; TBD
    'H2O':       unity_optical_depth = 9.e99 ; TBD
  ENDCASE
  
  ;===================================Solar Line Spectra===================================================
  CASE Line OF  
    'Na-D': begin
        skipline = 289106 ; 588.2 nm
        numline  = 2400   ; 590.6 nm
        spectrum = 'Kurucz_2005_irradthuwl.dat'   ; Great resolution. WL in nm, flux in flux is in W/m2/nm at 1AU
      end
    'K-D':  begin  
        skipline = 466901 ; 766.0 nm
        numline  = 4500   ; 770.5 nm
        spectrum = 'Kurucz_2005_irradthuwl.dat'   ; Great resolution. WL in nm, flux in flux is in W/m2/nm at 1AU
      end
    'Mg-2853': begin                              ; Need to use a lower resolution speectrum here. Should 
        spectrum = 'sao2010.solref.converted.txt' ; 0.4 A resolution (poor!), converted to 0.1A intervals. WL in nm, flux in Photons s-1 cm-2 nm-1 at 1AU
      end
  ENDCASE
  
  CASE spectrum OF  
    'Kurucz_2005_irradthuwl.dat': begin
        ; Read in a solar spectrum
        SOLAR_SPECTRUM_FILE = directory + '\Solar_Spectral_Irradiance\Kurucz\' + spectrum  ; Solar Spectrum at 1AU in W/m2/nm
        READCOL, SOLAR_SPECTRUM_FILE, F='F,F', WL_nm, flux, STRINGSKIP = '#', /NAN, /Silent, Numline = numline, skipline = skipline ;flux is in W/m2/nm
    
        ; change flux units from W/m^2/nm to photons / (cm^2 s A)
        ; multiply by ((lambda / hc) / 1 W)  * (1 m^2 / 1e4 cm^2) * (1 nm / 10 A)
        conversion = ((WL_nm*1.e-9)/(6.62606957e-34*299792458.D)) * (1./1.e4) * (1./10.)
        flux = flux * conversion                                                                    ; photons / (cm^2 s A)
        WL_A = temporary(WL_nm) * 10.                                                               ; Wavelength from nm into angstroms
      end
    'sao2010.solref.converted.txt': begin
        SOLAR_SPECTRUM_FILE = directory + '\Solar_Spectral_Irradiance\Kurucz\' + 'sao2010.solref.converted.txt'  ; WL in nm, flux in Photons s-1 cm-2 nm-1 at 1AU
        READCOL, SOLAR_SPECTRUM_FILE, F='F,F', WL_nm, flux, STRINGSKIP = '#', /NAN, /Silent 
        flux = flux / 10.                                                                           ; photons / (cm^2 s nm) to photons / (cm^2 s A)
        WL_A = temporary(WL_nm) * 10.                                                               ; Wavelength from nm into angstroms
      end
  ENDCASE
 
  Line_data = {name:name, line:line, unity_optical_depth:unity_optical_depth, wavelength:WL_A, intensity:flux}
  return
END