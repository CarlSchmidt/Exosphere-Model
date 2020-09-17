Pro load_particle_data, test_particle, line

COMMON Model_shared, Body, Ephemeris_time, Seed, Directory, Particle_data, Solar_data, Debug

;===================================FORMAT NAME STRINGS===============================================

CASE test_particle OF 
  'H': NAME = 'H'
  'O': NAME = 'O'
  'Na': NAME = 'Na'
  'Mg': NAME = 'Mg'
  'K': NAME = 'K'
  'H2O': NAME = 'H!D2!NO' 
ENDCASE

;===========================================PARTICLE MASSES======================================================= 

amu = 1.660538921e-27 ;one amu, atomic mass unit in kg
CASE test_particle OF 
  'H': mass = 1.00794*amu
  'O': mass = 15.9994*amu 
  'Na': mass = 22.98977*amu
  'Mg': mass = 24.305*amu
  'K': mass = 39.0983*amu
  'H2O': mass = 18.01528*amu
ENDCASE

;=========================================IONIZATION POTENTIALS===================================================

;cf. http://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
CASE test_particle OF ;units are eV 
  'H': Ionization_potential = 13.598433770784
  'O': Ionization_potential = 13.618054
  'Na': Ionization_potential = 5.1390766
  'Mg': Ionization_potential = 7.6462
  'K': Ionization_potential = 4.3407
  'H2O': Ionization_potential = 5.11359 ;H20 threshold is for all dissociation in to OH and H (Herzberg, 1966)
                                        ;12.65eV is H2O--> H2O+ (cf. http://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Units=SI&Mask=20)
ENDCASE

;=====================================PHOTO-IONIZATION CROSS-SECTIONS===============================================

CASE test_particle OF
  'H': begin
    phot_ioniz_X_sec = strcompress(directory + '\Cross-Sections\Photo\Huebner_2011\H_to_H+.txt') 
    if !VERSION.OS_FAMILY eq ('unix' or 'MacOS') then phot_ioniz_X_sec = repstr(phot_ioniz_X_sec , '\', '/')
    READCOL, phot_ioniz_X_sec, Format='A,A', Cross_wavelength, photo_cross_sec, /SILENT  
    start = max(where(Cross_wavelength eq 'Lambda'))+1
    Cross_wavelength = Cross_wavelength[start:*] / 10. ;convert from Angstoms to nm, solar irradiance from SEE/SORCE are in nm
    photo_cross_sec = photo_cross_sec[start:*] ;in cm^2
    if keyword_set(debug) then plot, Cross_wavelength, photo_cross_sec, /ylog, /xlog, Title = 'H Photo-Ionization Cross-Section',$
      xtitle = 'nm', ytitle = 'cm!u2!n'
  end
  'Na': begin
    phot_ioniz_X_sec = strcompress(directory + '\Cross-Sections\Photo\Huebner_2011\Na_to_Na+_theoretical.txt')
    if !VERSION.OS_FAMILY eq ('unix' or 'MacOS') then phot_ioniz_X_sec = repstr(phot_ioniz_X_sec , '\', '/')
    READCOL, phot_ioniz_X_sec, Format='A,A', Cross_wavelength, photo_cross_sec, /SILENT  
    start = max(where(Cross_wavelength eq 'Lambda'))+1
    Cross_wavelength = Cross_wavelength[start:*] / 10. ;convert from Angstoms to nm, solar irradiance from SEE/SORCE are in nm
    photo_cross_sec = photo_cross_sec[start:*] ;in cm^2
    if keyword_set(debug) then plot, Cross_wavelength, photo_cross_sec, /ylog, /xlog, Title = 'Theoretical Na Photo-Ionization Cross-Section',$
      xtitle = 'nm', ytitle = 'cm!u2!n'
    
    ;CS Note: Huebner's cross-sections when binned to his solar flux wavelengths seem high. 
    ;This suggests it's ok for my calculated photo-ionization lifefimes to be ~25% longer than Huebner et al., 1992
    ;I have checked Chang & Kelly (1975) for the original theoretical cross-section against his values in Huebner, 2011
    ;However, Chang and Kelly only quote cross-sections as a function of ejected electron momentum over a range .1-11.2 a.u.,
    ;simply equating this to incident photon momentum is not correct, i.e. lambda not equal (6.62606957e-34 / (((findgen(12)+1.)/10.) * 2.e-24) ) * 1.e10
    ;In short, the validity of my photo-ionization calculation hinges on Huebner's intepretation of Chang & Kelly (1975)  
    view_discrepancy = 1
    if view_discrepancy then begin
      phot_ioniz_X_sec = strcompress(directory + '\Cross-Sections\Photo\Huebner_2011\Na_to_Na+_theoretical.txt')
      if !VERSION.OS_FAMILY eq ('unix' or 'MacOS') then phot_ioniz_X_sec = repstr(phot_ioniz_X_sec , '\', '/')
      READCOL, phot_ioniz_X_sec, Format='A,A', Cross_wavelength, photo_cross_sec, /SILENT  
      start = max(where(Cross_wavelength eq 'Lambda'))+1
      Cross_wavelength = Cross_wavelength[start:*] / 10. ;convert from Angstoms to nm, solar irradiance from SEE/SORCE are in nm
      photo_cross_sec = photo_cross_sec[start:*] ;in cm^2
      if keyword_set(debug) then plot, Cross_wavelength, photo_cross_sec, Title = 'Theoretical Na Photo-Ionization Cross-Section',$
        xtitle = 'nm', ytitle = 'cm!u2!n', psym = 1., yrange = [5.e-21, 6.e-20], xrange = [200.,250.], Color=cgColor('black'), $
        Background=cgColor('white'), thick = 2., charthick = 1.6, charsize = 1.6
      ;Huebner's cross-sections when binned to his solar flux wavelengths seem high. 
      ;Plot his interpolated binned cross-sections in RED, and his unbinned data as quoted from Chang & Kelly (1975)
        phot_ioniz_X_sec = strcompress(directory + '\Cross-Sections\Photo\Huebner_2011\Na_to_Na+_theoretical_lambda_binned.txt')
        if !VERSION.OS_FAMILY eq ('unix' or 'MacOS') then phot_ioniz_X_sec = repstr(phot_ioniz_X_sec , '\', '/')
        READCOL, phot_ioniz_X_sec, Format='A,A', Cross_wavelength_binned, photo_cross_sec_binned, /SILENT  
        start = max(where(Cross_wavelength_binned eq 'Lambda'))+1
        Cross_wavelength_binned = Cross_wavelength_binned[start:*] / 10. ;convert from Angstoms to nm, solar irradiance from SEE/SORCE are in nm
        photo_cross_sec_binned = photo_cross_sec_binned[start:*] ;in cm^2
        if keyword_set(debug) then oplot, Cross_wavelength_binned, photo_cross_sec_binned , psym = 2., Color=cgColor('Red')        
     ;Indeed his interpolation seems inaccurate, but this still doesn't explain the full ~30% discrepancy from my value.
     ;At 1 AU Fulle et al. predicted 5.26e-6 /s, Huebner's best accepted value is 5.92e-6 /s, I get ~4.23e-6 /s  WHAT DOES THE MERCURY TAIL FALLOFF SHOW?
     endif
   end
   'Mg': begin    
     phot_ioniz_X_sec = strcompress(directory + '\Cross-Sections\Photo\Huebner_2011\Mg_to_Mg+.txt')
     if !VERSION.OS_FAMILY eq ('unix' or 'MacOS') then phot_ioniz_X_sec = repstr(phot_ioniz_X_sec , '\', '/')
     READCOL, phot_ioniz_X_sec, Format='A,A', Cross_wavelength, photo_cross_sec, /SILENT  
     start = max(where(Cross_wavelength eq 'Lambda'))+1
     Cross_wavelength = Cross_wavelength[start:*] / 10. ;convert from Angstoms to nm, solar irradiance from SEE/SORCE are in nm
     photo_cross_sec = photo_cross_sec[start:*] ;in cm^2
     if keyword_set(debug) then plot, Cross_wavelength, photo_cross_sec, /ylog, /xlog, Title = 'Mg Photo-Ionization Cross-Section',$
       xtitle = 'nm', ytitle = 'cm!u2!n'
   end  
   'K': begin
     phot_ioniz_X_sec = strcompress(directory + '\Cross-Sections\Photo\Huebner_2011\K_to_K+.txt')
     if !VERSION.OS_FAMILY eq ('unix' or 'MacOS') then phot_ioniz_X_sec = repstr(phot_ioniz_X_sec , '\', '/')
     READCOL, phot_ioniz_X_sec, Format='A,A', Cross_wavelength, photo_cross_sec, /SILENT     
     start = max(where(Cross_wavelength eq 'Lambda'))+1
     Cross_wavelength = Cross_wavelength[start:*] / 10. ;convert from Angstoms to nm, solar irradiance from SEE/SORCE are in nm
     photo_cross_sec = photo_cross_sec[start:*] ;in cm^2
     if keyword_set(debug) then plot, Cross_wavelength, photo_cross_sec, /ylog, /xlog, Title = 'Mg Photo-Ionization Cross-Section',$
       xtitle = 'nm', ytitle = 'cm!u2!n'
   end
   'H2O': begin
      phot_ioniz_X_sec = strcompress(directory + '\Cross-Sections\Photo\Huebner_2011\H2O_All.txt')
      if !VERSION.OS_FAMILY eq ('unix' or 'MacOS') then phot_ioniz_X_sec = repstr(phot_ioniz_X_sec , '\', '/')
      READCOL, phot_ioniz_X_sec, Format='A,A,A,A,A,A,A,A,A', Cross_wavelength, photo_cross_sec, H_OH, H2_O1D, O_H_H, OHplus_H, Oplus_H2, Hplus_OH, H2Oplus,  /SILENT  
      start = max(where(Cross_wavelength eq 'Lambda')) + 1
      Cross_wavelength = Cross_wavelength[start:*] / 10. ;convert from Angstoms to nm, solar irradiance from SEE/SORCE are in nm
      photo_cross_sec = photo_cross_sec[start:*] ;total combined cross-sections for all photo-decompositions in cm^2
      H_OH = H_OH[start:*] ;total combined cross-sections for all photo-decompositions in cm^2
      H2_O1D = H2_O1D[start:*] ;total combined cross-sections for all photo-decompositions in cm^2
      O_H_H = O_H_H[start:*] ;total combined cross-sections for all photo-decompositions in cm^2
      OHplus_H = OHplus_H[start:*] ;total combined cross-sections for all photo-decompositions in cm^2
      Oplus_H2 = Oplus_H2[start:*] ;total combined cross-sections for all photo-decompositions in cm^2
      Hplus_OH = Hplus_OH[start:*] ;total combined cross-sections for all photo-decompositions in cm^2
      H2Oplus = H2Oplus[start:*] ;total combined cross-sections for all photo-decompositions in cm^2
      ;photo_cross_sec = H_OH
      if keyword_set(debug) then plot, Cross_wavelength, photo_cross_sec, /ylog, /xlog, Title = 'H20 Photo-Ionization Cross-Section',$
        xtitle = 'nm', ytitle = 'cm!u2!n', charsize = 1.6 
      Result = DIALOG_MESSAGE('H2O ANTI-SUNWARD ACCELERATION IS TO ARBITRARILY SET 20 CM / S^2 at 1 AU')  
   end
endcase

photo_ionize_data = [transpose(Cross_wavelength), transpose(photo_cross_sec)]
Particle_data     = {Name:Name, Mass:Mass, Ionization_potential:Ionization_potential, photo_ionize_data:photo_ionize_data}

return
end