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
  COMMON Model_shared, Body, Ephemeris_time, Seed, Directory, Particle_data, Debug
  k = 1.3806504e-23      ;Boltzmann's constant in MKS units
  m = particle_data.mass ;Mass of the test particle in kg
  temperature = distribution_parameters.temperature ;Temperature passed from the main level input string
  f = 4.*!pi*(v^2.)*( m / (2.*!pi*k*temperature))^1.5 * exp((-m*v^2.)/(2.*k*temperature)) ;the distrubtion function over velocity, mks units
  f = f / max(f)
  return, f
end

function MBF, v, distribution_parameters
  COMMON Model_shared, Body, Ephemeris_time, Seed, Directory, Particle_data, Debug
  k = 1.3806504e-23      ;Boltzmann's constant in MKS units
  m = particle_data.mass ;Mass of the test particle in kg
  temperature = distribution_parameters.temperature ;Temperature passed from the main level input string
  f = !pi*(v^3.)*( m / (2.*!pi*k*temperature))^1.5 * exp((-m*v^2.)/(2.*k*temperature)) ;the distrubtion function over velocity, mks units
  ;pi factor: integral of sin(theta)cos(theta)dtheta from 0 to pi/2 (outward flow only) is .5, and integral of dphi from 0 to 2pi is 2Pi. But it of course doesn't matter
  f = f / max(f)
  return, f
end

function Kappa, v, distribution_parameters
  COMMON Model_shared, Body, Ephemeris_time, Seed, Directory, Particle_data, Debug
  k = 1.3806504e-23      ;Boltzmann's constant in MKS units
  m = particle_data.mass ;Mass of the test particle in kg
  temperature = distribution_parameters.temperature ;Temperature passed from the main level input string
  kappa = distribution_parameters.kappa
  
  ;the following is from Summers and Thorne (1991) and applies to velocity, not speed as it should and it differs from Moore and Mendillo, 2005
  ;theta = sqrt(((2.*kappa)-3.)/kappa)*sqrt(k*temperature/m)
  ;distfunct_v = (Normal/(theta*sqrt(!pi))) * (gamma(kappa+1.)/((kappa^1.5)*gamma(kappa-.5))) * (1.+((v_array^2.)/(kappa*(theta)^2.)))^(-kappa)
  
  ;from Moore and Mendillo, 2005:
  ;distfunct_v =  4.*!pi*(v_array^2.) * (Kappa*m/(((2.*kappa)-3.)*!pi*k*Temperature))^1.5 * $
  ;              (gamma(kappa+1.)/((kappa^1.5)*gamma(kappa-.5))) * (1. + (m*kappa*(v_array^2.))/(((2.*kappa)-3.)*k*Temperature))^(-(kappa+1.))
  
  ;The review by Pierrard and Lazar, 2010 differs from both of the above:
  w = sqrt(((2.*kappa)-3.)*k*temperature/(kappa*m)) 
  f =  4.*!pi*(v^2.) * 1./(2.*!pi*(kappa*w^2.)^1.5) * (gamma(kappa+1.)/(gamma(kappa-.5)*gamma(1.5))) * (1. + (v^2.)/(kappa*w^2.))^(-(kappa+1.))
  ;f is the distrubtion function over velocity, mks units
  
  ;Pierrard and Lazar, 2010 is the one that best matches the corresponding Maxwell-Bolzmann distribution
  f = f / max(f)
  return, f
end

function Shematovich, v, distribution_parameters ;
  COMMON Model_shared, Body, Ephemeris_time, Seed, Directory, Particle_data, Debug
  PDF = distribution_parameters.PDF
  Velocity = distribution_parameters.velocity
  f = interpol(PDF, Velocity, v, /NaN)
  f = f / max(f)
  return, f
end

function generate_velocity_distribution, number_of_particles, distribution_parameters, Force_distribution = Force_distribution  
  
  COMMON Model_shared, Body, Ephemeris_time, Seed, Directory, Particle_data, Debug
  
  k = 1.3806504e-23      ;Boltzmann's constant in MKS units
  m = particle_data.mass ;Mass of the test particle in kg
  distribution = distribution_parameters.distribution ;functional form of the probability distribution
  
  if not keyword_set(Force_distribution) then begin
      temperature = distribution_parameters.temperature
      V_therm = sqrt(2.*k*temperature/m)  ; thermal velocity    
      vmax = 4.*V_therm  ;in a practical Monte Carlo sense, the maximum velocity possible for a given temperature
      xrange = [0.,vmax]
      yrange = [0.,1.]
      delx = max(xrange)-min(xrange) 
      ;;; Set up some useful arrays
      bad = lindgen(number_of_particles)
      pxi = fltarr(number_of_particles)
      xi = fltarr(number_of_particles)
      yi = fltarr(number_of_particles)
      nbad = number_of_particles
      cool = 0b                       ;loop control variable
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
      ; -> randomize them so that they differ for each particle follwing the random number seed.
      randomNumbers = RANDOMU(seed, number_of_particles) 
      scrambled_indicies = SORT(randomNumbers)
      speed = xi[scrambled_indicies]
      v_array = findgen(number_of_particles)/float(number_of_particles)*vmax
      distfunct_v = call_function(distribution, v_array, distribution_parameters)
  endif
  
  if keyword_set(Force_distribution) then begin
    ;******************Maxwell-Boltzmann and Generalized Lorentzian (Kappa) Distributions*********************************************************
    if ((STRPOS(distribution, 'Maxwellian') eq 0) or (STRPOS(distribution, 'Kappa') eq 0)) then begin
      temperature = distribution_parameters.temperature
      V_therm = sqrt(2.*k*temperature/m)  ; thermal velocity    
      vmax = 4.*V_therm  ;in a practical Monte Carlo sense, the maximum velocity possible for a given temperature
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
    cghistoplot, speed/1000., Xtitle = 'Particle Release Velocity Distribution (km/s)', ytitle = $ 
      strcompress('Number per'+string(number_of_particles)+' Test Particles'), $
      charsize=1.6, histdata = histdata, binsize = .1
      
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

return, speed/1000. ;return speeds in km/s rather then m/s
end