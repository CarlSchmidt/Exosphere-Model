; NAME:
;        generate_velocity_distribution
;
; PURPOSE:
;        Draw a random set of N values with velocity V following a probability
;        distribution function P(V).
;
; CALLING SEQUENCE:
;
;        result = generate_velocity_distribution(number_of_particles, distribution_parameters [, Force_distribution])
;
; INPUTS:
;
;        number_of_particles      -  The number of random values to be returned
;
;        distribution_parameters  -  An stucture containing the parameters of the 
;                                    associated distribution function    
    
; OUTPUTS:
;
;        result                   -  A random set of number_of_particles values of V which follow
;                                    the PDF defined in distribution.
;
; KEYWORDS:
; 
;        Force_distribution       -  The generates a perfect distribution, but one
;                                    with non-unique velocities. Output has
;                                    randomized ordering. The order will differ as
;                                    the seed updates, but the velocities are not
;                                    random. 
;                                
;                                    The default is Monte Carlo velocities that have 
;                                    statistical deviations from a perfect distribution.
;                                    Such deviations approach zero for a large
;                                    number_of_particles, or over many iterations
;                                    provided seed is being correctly passed.
;
; PROCEDURE:
;           For a single value (number_of_particles = 1):
;
;             1) Choose a uniform random number Xi between XRANGE[0] and
;                XRANGE[1] and a uniform random value Yi between 0 and 1
;             2) Evaluate P(Xi)
;                  a) if P(Xi) > Yi then return Xi
;                  b) else return to step 1
;                  
; WRITTEN: 
;         Carl Schmidt, UVa 2013.           
;         CS 2021 Added: Shematovich (2013) velocity distribution

function Maxwellian, v, distribution_parameters
  COMMON Model_shared, Body, Ephemeris_time, Obs_Body_Ltime, Parent_ID, Seed, Directory, Particle_data, Line_data, Debug
  k = 1.3806504e-23      ;Boltzmann's constant in MKS units
  m = particle_data.mass ;Mass of the test particle in kg
  temperature = distribution_parameters.temperature ;Temperature passed from the main level input string
  f = 4.*!pi*(v^2.)*( m / (2.*!pi*k*temperature))^1.5 * exp((-m*v^2.)/(2.*k*temperature)) ;the distrubtion function over velocity, mks units
  f = f / max(f)
  return, f
end

function Johnson_2002, v, distribution_parameters
  COMMON Model_shared, Body, Ephemeris_time, Obs_Body_Ltime, Parent_ID, Seed, Directory, Particle_data, Line_data, Debug
  E = 0.5*Particle_data.mass*(v^2)     ; Joules
  E = E * 6.242e18                     ; convert Joules to eV 
  x = distribution_parameters.x_param
  U = distribution_parameters.U_param  ; U, units of eV
  f = E * (U^x) /( (E + U)^(2. + x) )
  f = f / max(f)
  return, f
end
 
function MBF, v, distribution_parameters
  COMMON Model_shared, Body, Ephemeris_time, Obs_Body_Ltime, Parent_ID, Seed, Directory, Particle_data, Line_data, Debug  
  k = 1.3806504e-23                    ; Boltzmann's constant in MKS units
  m = particle_data.mass               ; Mass of the test particle in kg
  temperature = distribution_parameters.temperature ;Temperature passed from the main level input string
  f = !pi*(v^3.)*( m / (2.*!pi*k*temperature))^1.5 * exp((-m*v^2.)/(2.*k*temperature)) ; the distrubtion function over velocity, mks units
                                                                                       ; pi factor in front is for completeness here:
                                                                                       ; integral of sin(theta)cos(theta)dtheta from 0 to pi/2 (outward flow only) is .5
                                                                                       ; integral of dphi from 0 to 2pi is 2Pi. 
                                                                                       ; Of course this pi factor does doesn't matter since it's normalized.
  f = f / max(f)
  return, f
end

function Kappa, v, distribution_parameters
  COMMON Model_shared, Body, Ephemeris_time, Obs_Body_Ltime, Parent_ID, Seed, Directory, Particle_data, Line_data, Debug
  k = 1.3806504e-23      ;Boltzmann's constant in MKS units
  m = particle_data.mass ;Mass of the test particle in kg
  temperature = distribution_parameters.temperature ;Temperature passed from the main level input string
  kappa = distribution_parameters.kappa
  
  ; Note that there is some ambiguity about what people mean when referencing to a Kappa distribution.
  
  ; The following is from Summers and Thorne (1991) and applies to velocity, not speed as it should and it differs from Moore and Mendillo, 2005
  ; theta = sqrt(((2.*kappa)-3.)/kappa)*sqrt(k*temperature/m)
  ; distfunct_v = (Normal/(theta*sqrt(!pi))) * (gamma(kappa+1.)/((kappa^1.5)*gamma(kappa-.5))) * (1.+((v_array^2.)/(kappa*(theta)^2.)))^(-kappa)
  
  ; From Moore and Mendillo, 2005:
  ; distfunct_v =  4.*!pi*(v_array^2.) * (Kappa*m/(((2.*kappa)-3.)*!pi*k*Temperature))^1.5 * $
  ;              (gamma(kappa+1.)/((kappa^1.5)*gamma(kappa-.5))) * (1. + (m*kappa*(v_array^2.))/(((2.*kappa)-3.)*k*Temperature))^(-(kappa+1.))
  
  ; The review by Pierrard and Lazar, 2010 differs from both of the above:
    w = sqrt(((2.*kappa)-3.)*k*temperature/(kappa*m)) 
    f =  4.*!pi*(v^2.) * 1./(2.*!pi*(kappa*w^2.)^1.5) * (gamma(kappa+1.)/(gamma(kappa-.5)*gamma(1.5))) * (1. + (v^2.)/(kappa*w^2.))^(-(kappa+1.))
    ;f is the distrubtion function over velocity, mks units
  
  ; Pierrard and Lazar, 2010 is the one that best matches the corresponding Maxwell-Bolzmann distribution
    f = f / max(f)
    return, f
end

function Shematovich, v, distribution_parameters ;
  COMMON Model_shared, Body, Ephemeris_time, Obs_Body_Ltime, Parent_ID, Seed, Directory, Particle_data, Line_data, Debug
  PDF      = distribution_parameters.PDF
  Velocity = distribution_parameters.velocity
  f        = interpol(PDF, Velocity, v, /NaN)
  f        = f / max(f)
  return, f
end

function Sputtering_ST, v, distribution_parameters ; Sigmund Thompson,
  ; See Wurz et al. 2010 for Equations applied 
  ; See Morrissey et al. 2022 Fig. 3 for sodium energy distribution 
  COMMON Model_shared, Body, Ephemeris_time, Obs_Body_Ltime, Parent_ID, Seed, Directory, Particle_data, Line_data, Debug
  SBE = distribution_parameters.SBE
  m   = particle_data.mass                        ;  Mass of the test particle in kg
  E_e = (1./2.)*m*(v^2.) / 1.60217634e-19                 ; Indep variable: Energy converted from Joules to eV
  E_i = 1.e3                                              ; Energy of the Projectile particle in eV, use 1 keV protons for an estimate
  E_b = SBE                                               ; Surface binding energy in units of eV
  m = m*1.e10                                             ; Multiply mass x 10^10 to Avoid NAN in E_c, this arbitrary factor cancels out in E_c
  M_projectile = 1.672621777e-27 * 1.e10                  ; Mass of a proton in kg (x 10^10 to Avoid NAN in E_c) use 1kev protons for an estimate
  E_c = E_i*(4.*M*M_projectile)/((M+M_projectile)^2.)     ; Cutoff energy: the maximum energy a projectile can give to the atom
  f = (6.*E_b/(3.-(8.*sqrt(E_b/E_c)))) * E_e/((E_e+E_b)^3.) * (1.-sqrt((E_e+E_b)/E_c)) ;Energy distribution for Ion Sputtering (Wurz et al., 2010)
  f = f / max(f)
;  cgplot, e_e, f, psym=4, /xlog, xr = [0.1,40], /xs, $   ; Used to Verify that e vs f mataches Morrissey et al. 2022 Fig. 3 
;    xtit = 'Na atom ejection energy [eV]', ytit = 'Energy Distribution Function [dN/N]', charsize = 3
  return, f
end

function generate_velocity_distribution, number_of_particles, distribution_parameters, Force_distribution = Force_distribution  
  COMMON Model_shared, Body, Ephemeris_time, Obs_Body_Ltime, Parent_ID, Seed, Directory, Particle_data, Line_data, Debug
  
  k = 1.3806504e-23      ; Boltzmann's constant in MKS units
  m = particle_data.mass ; Mass of the test particle in kg
  distribution = distribution_parameters.distribution ;functional form of the probability distribution
  
  ; First, we'll need to establish the span the VDF range, used for xrange below
  ; The span should just cover the VDF for a practial number of test particles (say 100,000)
  ; Too wide a range works, but is computationally inefficient. 
  ; So in a practical Monte Carlo sense, get a maximum velocity possible for the 0.00001 upper range of the specified VDF 
  
  Case distribution of 
    'Maxwellian'   : begin
                      temperature = distribution_parameters.temperature
                      V_therm = sqrt(2.*k*temperature/m)  ; thermal velocity    
                      vmax = 4.*V_therm                   ; distribution apporaches 0.00001 here
                    end
    'MBF'          : begin
                      temperature = distribution_parameters.temperature
                      V_therm = sqrt(2.*k*temperature/m)  ; thermal velocity    
                      vmax = 4.5*V_therm                  ; distribution apporaches 0.00001 here
                    end
    'Kappa'        : begin
                      temperature = distribution_parameters.temperature
                      V_therm = sqrt(2.*k*temperature/m)  ; thermal velocity
                      vmax = 5.0*V_therm                  ; distribution apporaches 0.00001 here
                    end                
    'Johnson_2002' : Vmax = 10000.                        ; 10 km/s, corresponds to ~6eV for Na, distribution approaches 0.00001 here
    'Sputtering_ST': Vmax = 32767.                        ; 32 km/s, corresponds to ~128eV for Na, sputtering distribution approaches 0.00084 here
  endcase
  
  if not keyword_set(Force_distribution) then begin
      xrange = [0.,vmax]
      yrange = [0.,1.]
      ;;; Set up some useful arrays
      bad = lindgen(number_of_particles)
      pxi = fltarr(number_of_particles)
      xi = fltarr(number_of_particles)
      yi = fltarr(number_of_particles)
      nbad = number_of_particles
      cool = 0b                           ; loop control variable
      ct = 0
      while not cool do begin
          ;;; Get set of x values for values not yet settled
          xi[bad] = randomu(seed,nbad)*(max(xrange)-min(xrange))+min(xrange)
          ;;; Get set of y values for values not yet settled
          yi[bad] = randomu(seed,nbad)*(max(yrange)-min(yrange))+min(yrange)
          ;;; Get set of comparison values for indices not yet settle
          pxi[bad] = call_function(distribution, xi[bad], distribution_parameters)
          ;;; Is pxi > yi?  If not, mark the elements which require
          ;;; changing and rinse, repeat. Dry once there are no "bad"
          ;;; values 
          bad = where(yi gt pxi, nbad, complement = good)
          ngood = n_elements(good)
          ct += 1.
          if nbad eq 0 then cool = 1
      endwhile
      
      ;speed now has all the speeds that will be assigned to the particles.
      ; -> randomize them so that they differ for each particle following the random number seed.
      randomNumbers      = RANDOMU(seed, number_of_particles) 
      scrambled_indicies = SORT(randomNumbers)
      speed              = xi[scrambled_indicies]
      v_array            = findgen(number_of_particles)/float(number_of_particles)*vmax
      distfunct_v = call_function(distribution, v_array, distribution_parameters)
  endif
  
  if keyword_set(Force_distribution) then begin ; not sure all this is well tested.
    ;******************Maxwell-Boltzmann and Generalized Lorentzian (Kappa) Distributions*********************************************************
    if ((STRPOS(distribution, 'Maxwellian') eq 0) or (STRPOS(distribution, 'Kappa') eq 0)) then begin
      temperature = distribution_parameters.temperature
      V_therm = sqrt(2.*k*temperature/m)  ; thermal velocity    
      vmax = 4.*V_therm                   ; in a practical Monte Carlo sense, the maximum velocity possible for a given temperature
      v_array = findgen(number_of_particles)/float(number_of_particles)*vmax
      distfunct_v = call_function(distribution, v_array, distribution_parameters)
    endif
    ;*****************************************************************************************************************************************
    if keyword_set(distfunct_v) then begin ;for distribution functions with WRT speed
      ;if the probability distribution is normalized, it integrates to one.
      sumdist_v = distfunct_v ;sum up the distribtion function element by element into an array
      for i = 1, n_elements(v_array)-2 do begin
          sumdist_v(i) = sumdist_v(i-1)+((distfunct_v(i-1) + distfunct_v(i))/2.)
      endfor
      sumdist_v = sumdist_v/max(sumdist_v) ;normalize to a maximum of one
      distfunct_v = distfunct_v /float(max(sumdist_v))
      r = findgen(number_of_particles)/number_of_particles ; an number_of_particles length array of increasingly large fractions of unity
      speed = fltarr(number_of_particles) ;speed is the speed each particle will have 
      i = long(0)
      while i le n_elements(sumdist_v)-2 do begin ;step through the sum of the distribution 
          top = sumdist_v(i)      
          bottom = 0.
          if i gt 0 then bottom = sumdist_v(i-1) ;the lower bound of the sum of the distribution 
          inrange = where((r lt top) and (r ge bottom))
          if inrange(0) ne -1 then speed[inrange] = i
          i = i + 1L
          if sumdist_v(i-1) eq 1. then i = n_elements(sumdist_v )- 1 ; when stepping through the distribution totals to 1, stop the loop
      endwhile  
      speed = speed*(VMAX/FLOAT(number_of_particles)) ;re-scale to the correct speed
    
      ;speed now has all the speeds that will be assigned to the particles.
      ; -> randomize them so that they differ for each particle follwing the random number seed.
      randomNumbers = RANDOMU(seed, number_of_particles) 
      scrambled_indicies = SORT(randomNumbers)
      speed[indgen(number_of_particles)] = speed[scrambled_indicies]
    endif  
  endif
  ;************************************************************************************************************************************************************
  if keyword_set(debug) then begin 
    window, 1
    cghistoplot, speed/1000., Xtitle = 'Particle Release Velocity Distribution (km/s)', ytitle = $
      strcompress('Number per'+string(number_of_particles)+' Test Particles'), $
      charsize=1.6, histdata = histdata, binsize = .1

    ; Overplot the velocity distribution function to inspect the match betweeen the Monte Carlo histogram and the analytic function
      overplot_velocities      = findgen(fix(vmax)) ; in m/s
      
      Case distribution of
        'Maxwellian'   : distribution_to_overplot = Maxwellian(overplot_velocities, distribution_parameters)
        'MBF'          : distribution_to_overplot = MBF(overplot_velocities, distribution_parameters)
        'Johnson_2002' : distribution_to_overplot = Johnson_2002(overplot_velocities, distribution_parameters)
        'Sputtering_ST': distribution_to_overplot = Sputtering_ST(overplot_velocities, distribution_parameters)
      endcase
      cgplot, overplot_velocities/1000., distribution_to_overplot*max(histdata), /overplot ; plot in km/s units, normalize the analytic function to the histogram peak 
      
    ; If we're simulating superthermal Hydrogen at Mars, inspect the VDF specific to that population   
      if distribution eq 'Shematovich' then begin
        readcol, directory+'Shematovich_2013_Fig5.txt', F='A,A', v1, v2, DELIMITER = ','
        energy = float(v1) ;eV
        PDF = float(v2) ;cm-2, s-1, eV-1
        Velocity = sqrt(2.*energy*1.60218e-19/1.67377e-27) ;convert hydrogen energy in eV to joules to velocity in m/s
        eV_over_m_per_s = 1. / sqrt(2.*1.*1.60218e-19/1.67377e-27)
        PDF = PDF*eV_over_m_per_s ;convert hydrogen energy in eV to joules to velocity in m/s
        cgplot, Velocity/1000., PDF*max(histdata)/max(pdf), xtitle = '(H atom velocity km/s)', ytitle = 'cm-2, s-1, (m/s)-1 at 400 km', /overplot
      endif
  endif

return, speed/1000.  ; returns the distribution of particle speeds in km/s
end