PRO generate_line_profiles, Model_Image_RV, Model_Image_T=Model_Image_T

  ; REVISION HISTORY:
  ;    Written:                     P. Lierle, Jun 2026

  COMMON Model_shared, Body, Ephemeris_time, Obs_Body_Ltime, Parent_ID, Seed, Directory, Particle_data, Line_data, Debug
  COMMON Output_shared, Plot_range, Output_Size_In_Pixels, Output_Title, Center_in_frame, viewpoint, FOV, N_ticks, Ecliptic_Tickstep, Sun_Body_Tickstep, Observatory, Above_Ecliptic, $
    Sun_Body, Boresight_Pixel, Aperture_Corners, Measure_Linewidths, Inst_Pointing, Inst_Diameter, Inst_LSF_wav, Line
 
  k = 1.38064852D-16                                    ; cgs units erg * K^-1

  Morton_table_Na = $ ; See Morton et al., 2003 Table 4
    ;    Air       Vacuum   E_low cm-1  E_up cm-1   gl gu fl fu     f     loglamf error
    [[5895.9333D, 5897.5672, +0.022161, 16956.16631, 2, 2, 2, 1, 1.600D-01, 2.975, 33E-5],$
    [5895.9311D, 5897.5650, +0.022161, 16956.17261, 2, 2, 2, 2, 1.600D-01, 2.975, 33E-5],$
    [5895.9128D, 5897.5467, -0.036934, 16956.16631, 2, 2, 1, 1, 5.335D-02, 2.498, 33E-5],$
    [5895.9106D, 5897.5445, -0.036934, 16956.17261, 2, 2, 1, 2, 2.668D-01, 3.197, 33E-5],$

    [5889.9592D, 5891.5915, +0.022161, 16973.36448, 2, 4, 2, 1, 3.204D-02, 2.276, 33E-5],$
    [5889.9588D, 5891.5911, +0.022161, 16973.36563, 2, 4, 2, 2, 1.602D-01, 2.975, 33E-5],$
    [5889.9582D, 5891.5905, +0.022161, 16973.36757, 2, 4, 2, 3, 4.486D-01, 3.422, 33E-5],$
    [5889.9389D, 5891.5712, -0.036934, 16973.36396, 2, 4, 1, 0, 1.068D-01, 2.799, 33E-5],$
    [5889.9387D, 5891.5710, -0.036934, 16973.36448, 2, 4, 1, 1, 2.670D-01, 3.197, 33E-5],$
    [5889.9383D, 5891.5706, -0.036934, 16973.36563, 2, 4, 1, 2, 2.670D-01, 3.197, 33E-5]]

  case Line of
    'Na-D': begin
      m         = 22.989769D * 1.66053892D-24   ; cgs unit g
      ; Weight oscillator strengths by energy level population
      D1_f      = reform(Morton_table_Na[8,0:3])
      D1_f[0:1] = D1_f[0:1] * 5./8
      D1_f[2:3] = D1_f[2:3] * 3./8
      D2_f      = reform(Morton_table_Na[8,4:9])
      D2_f[0:2] = D2_f[0:2] * 5./8
      D2_f[3:5] = D2_f[3:5] * 3./8

      D1_struct = {name:'NaD1',               $
        lambda_air: 5895.9242D,    $
        f: D1_f,                   $
        lambda_0: reform(Morton_table_Na[0,0:3])}

      D2_struct = {name:'NaD2',               $
        lambda_air: 5889.9510D,    $
        f: D2_f,                   $
        lambda_0: reform(Morton_table_Na[0,4:9])}

      ; Structures allow for a different number of components per line
      line_info = {D1:D1_struct, D2:D2_struct}
    end
  endcase
  
  ; --------------------------------------------------- GENERATE EMISSION LINE PROFILES WITHIN INSTRUMENT APERTURE --------------------------------------------------------------------
  restore, strcompress(directory+Output_title+'Radial_v_Hist.sav')
  
  ; Plot radial velocity histogram in phot/s
  cgPS_Open, filename=strcompress(directory + Output_title + '_Radial_v_Hist.eps'), /ENCAPSULATED, /NOMATCH
    !P.font=1
    !p.charsize=1.5
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    cgplot, Velocity_Centers, Radial_v_Hist_avg, psym=10, xr=Velocity_Centers[minmax(where(Radial_v_Hist_avg gt 0.))], xtit='Radial Velocity [km/s]', ytit='Flux [photon/s]', $
      tit='Atom Radial Velocity Histogram', aspect=1.
  cgPS_Close

  ; Initialize arrays on which to build line profiles
  nlines          = n_elements(tag_names(line_info))
  lambda_binsize  = 0.001
  lambda          = fltarr(nlines, ceil(4/lambda_binsize))
  Line_Profile    = dblarr(nlines, ceil(4/lambda_binsize))
  for i=0, nlines-1 do begin
    transition    = line_info.(i)
    lambda[i,*]   = transition.lambda_air - 2.0 + lambda_binsize*indgen(ceil(4/lambda_binsize))
  endfor

  ; Instrumental LSF convolution can be applied either to the velocity histogram or the constructed line profile in
  ; wavelength space. However, instrumental LSF for a grating is constant in velocity space, so this is preferred for lines
  ; contained within a single order but spaced apart (i.e. sodium).
  Inst_LSF_vel  = cspice_clight() * (Inst_LSF_wav / transition.lambda_air)                ; Inst LSF sigma in km/s
  Inst_LSF_px   = Inst_LSF_vel / (Velocity_Centers[1]-Velocity_Centers[0])                ; Inst LSF sigma in velocity histogram pixels
  Radial_v_Hist_conv = convol(Radial_v_Hist_avg, gaussian_function(Inst_LSF_px), /edge_zero, /normalize)

  if total(Radial_v_Hist_conv) eq 0. then print, 'No particles within instrument aperture'

  ; Construct Line Profile
  for i=0, nlines-1 do begin
    transition  = line_info.(i)
    lambda_0    = transition.lambda_0
    f           = transition.f

    ; Add up fine/hyperfine components
    for j=0, lambda_0.length -1 do begin
      ; Determine for a hyperfine component the velocity that would dopplershift its center to the current pixel
      Doppler_v = cspice_clight() * (lambda[i,*]/lambda_0[j] - 1.)
      ; Interpolate the velocity histogram to get values at these doppler shift velocities
      v_at_lambda = interpol(Radial_v_Hist_conv, Velocity_Centers, Doppler_v)
      ; Construct profile for this hyperfine component including relative strength
      Line_Profile[i,*] += v_at_lambda * f[j]
    endfor

    ; Plot profile
    cgPS_Open, filename=strcompress(directory + Output_title + '_Line_Profile_'+transition.name+'.eps'), /ENCAPSULATED, xs=8, ys=6
    !P.font=1
    !p.charsize=1.5
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    cgplot, Lambda[i,*], Line_Profile[i,*], xr=Lambda[i,minmax(where(Line_Profile[i,*] gt 0.))], thick=3, xtit='Wavelength [A]', ytit='Flux [photon/s]', tit=transition.name+' Emission Line Profile', aspect=1.
    cgPS_Close
  endfor

  ; --------------------------------------------------- GENERATE EMISSION TEMPERATURE IMAGE FROM LINE PROFILES --------------------------------------------------------------------
  stackedModel_Image_CD = MRDFITS(strcompress(directory + Output_Title + '.fit'), 1, /SILENT)
  Model_Image_T         = fltarr(Output_Size_In_Pixels[0],Output_Size_In_Pixels[1])
  Model_Image_Terr      = fltarr(Output_Size_In_Pixels[0],Output_Size_In_Pixels[1])
  Model_Image_RV       /= mean(Model_Image_RV)    ; get numbers more reasonable, calibration doesn't matter

  for x=0, Output_Size_In_Pixels[0] -1 do begin
    for y=0, Output_Size_In_Pixels[1] -1 do begin

      ; Only include pixels with at least 10^8 atoms
      if stackedModel_Image_CD[x,y] lt 1.e8 then continue

      radial_v_hist = Model_Image_RV[x,y,*]
      T             = fltarr(nlines)
      Terr          = fltarr(nlines)

      ; Initialize arrays on which to build line profiles
      nlines          = n_elements(tag_names(line_info))
      lambda_binsize  = 0.001
      lambda          = fltarr(nlines, ceil(4/lambda_binsize))
      Line_Profile    = dblarr(nlines, ceil(4/lambda_binsize))
      for i=0, nlines-1 do begin
        transition    = line_info.(i)
        lambda[i,*]   = transition.lambda_air - 2.0 + lambda_binsize*indgen(ceil(4/lambda_binsize))
      endfor

      ; Construct Line Profile
      for i=0, nlines-1 do begin
        transition  = line_info.(i)
        lambda_0    = transition.lambda_0
        f           = transition.f

        ; Add up fine/hyperfine components
        for j=0, lambda_0.length -1 do begin
          ; Determine for a hyperfine component the velocity that would dopplershift its center to the current pixel
          Doppler_v = cspice_clight() * (lambda[i,*]/lambda_0[j] - 1.)
          ; Interpolate the velocity histogram to get values at these doppler shift velocities
          v_at_lambda = interpol(radial_v_hist, Velocity_Centers, Doppler_v)
          ; Construct profile for this hyperfine component including relative strength
          Line_Profile[i,*] += v_at_lambda * f[j]
        endfor

        ; Set up fitting parameters
        Line_Profile[i, where(Line_Profile[i,*] eq 0)] = !values.F_NaN    ; mask out everywhere but line
        max = max(Line_Profile[i,*], ind, /NaN)
        estimates = [lambda[i,ind], 0.02, max]
        parinfo               = replicate( {value:0., fixed: 0b, limited: [0b,0b], limits: dblarr(2), step:0}, 3)
        parinfo[0].limited    = 1b
        parinfo[0].limits     = estimates[0] + [-2., 2]               ; line center within 2A of estimate
        parinfo[1].limited[1] = 1b
        parinfo[1].limits[1]  = 1.0                                   ; sigma less than 1A
        parinfo[2].limited[0] = 1b
        parinfo[2].limits[0]  = 1.e-25                                ; amplitude positive
        parinfo[*].value      = [estimates[0],estimates[1],estimates[2]]

        ; Fit line
        sy = sqrt(Line_Profile[i,*])  ; poisson stats
        p = mpfitfun('GAUSS1', lambda[i,*], Line_Profile[i,*], sy, estimates, parinfo=parinfo, functargs={peak:1}, STATUS=status, yfit=fit, PERROR=perror, /NaN, /QUIET)
        sigma = p[1]

        ; Convert gaussian sigma to temperature
        Na_Temperature_Lookup, sigma, line=strmid(transition.name,2,4), error=perror[1], /sigma, temperature=temperature, errout=errout
        T[i]    = temperature
        Terr[i] = errout
      endfor

      ; Average the two lines
      Model_Image_T[x,y]    = avg(T)
      Model_Image_Terr[x,y] = avg(Terr)
    endfor
  endfor

  Model_Image_T_OG = Model_Image_T

  ; cut out temps below 0 as erroneous
  Model_Image_T[where(Model_Image_T le 0)] = !values.F_NaN

  ; cut out the top 1% of temps as erroneous
  sorted_T = Model_Image_T[where(finite(Model_Image_T))]
  sorted_T = sorted_T[sort(sorted_T)]
  threshold_T = sorted_T[sorted_T.length * 0.99]
  Model_Image_T[where(Model_Image_T gt threshold_T)] = !values.F_NaN
  
  RETURN
  
END