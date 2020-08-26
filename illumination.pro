Function Vector_Sep_angle, V1, V2 ;finds the angular separation between V2 & V1 (where unlike IDL V1 and V2 can be arrays of vectors!) 
  MAGNITUDE_v1 = SQRT( TOTAL( v1 * v1, 1 ) )
  MAGNITUDE_v2 = SQRT( TOTAL( v2 * v2, 1 ) )
  DOT_PRODUCT  = TOTAL( v1 * v2, 1 )
  return, ACOS( DOT_PRODUCT / (MAGNITUDE_v1*MAGNITUDE_v2) )
end

Function Vector_Sep_angle2, V1, V2 ;Alternate way to find the angular separation between V2 & V1 the way cspice_sep would do it 
  v3 = v1 - v2
  length = v3 / SQRT( TOTAL( V3 * V3, 1 ) )
  return, 2.d * asin( length/2.d )
end

function illumination, time, x, y, z

  ; Calculates the illumination on a bunch of points in space
  ; x,y,z: input location of points in km, relative to the body at the time (ephemeris time - time) in J2000 coordinates 
  ; Time is the number of seconds before the ephemeris time (Ephemeris_time) when the image is taken. 
  
  ; Assumes all bodies are circular and uses equatorial radius.
  ; Ignores Solar limb darkening
  
  COMMON Model_shared, Body, Ephemeris_time, Seed, Directory, Particle_data, Line_data, Debug
  
  cspice_bodvrd, 'Sun', 'RADII', 3, Solar_radius ; Find the planet's radius in Km
  cspice_bodvrd, Body , 'RADII', 3, Body_radius  ; Find the Body's radius in Km
  solar_radius = Solar_radius[0]                 ; Default to *Equatorial* radius of Sun in km
  Body_radius  = Body_radius[0]                  ; Default to *Equatorial* radius of the body in km

; Also Check for Potential Parent Body Occultations (Moons only)
  cspice_bods2c, body, bodyID_code, found
  radius_found = cspice_bodfnd( bodyID_code, "RADII" )
  if ((STRPOS(string(bodyID_code), '99') eq -1) and (strlen(strcompress(bodyID_code, /remove_all)) eq 3)) then begin 
    Parent_ID = strcompress(strmid(strcompress(bodyID_code, /remove_all),0,1) + '99')
    cspice_bodvrd, Parent_ID, 'RADII', 3, Parent_Radius ;Find G*Mass in in Km^3/s^2 
    Parent_radius = Parent_radius[0]
  endif

; Find the position of the body WRT the Sun for each particles time in units of km
  CSPICE_SPKPOS, body, ephemeris_time - REFORM(time), 'J2000', 'NONE', 'Sun', Sun2body, ltime  ; HACK----Not sure if light time is importnant here!!!
  sun_particle        = Sun2body[0:2,*] + [x,y,z]  ;Calculate the vectors from the Sun to the particles
  dist2body           = SQRT( TOTAL( [x,y,z] * [x,y,z], 1 ) )
  dist2Sun            = SQRT( TOTAL( sun_particle * sun_particle, 1 ) )
  angular_radius_body = asin(body_radius / dist2body) 
  angular_radius_Sun  = asin(solar_radius / dist2Sun) 

  ; If the body is a moon, then also find the position of Parent WRT the Sun for each particles time in units of km
    if keyword_set(Parent_ID) then begin
      CSPICE_SPKPOS, Parent_ID, ephemeris_time - REFORM(time), 'J2000', 'NONE', 'Sun', Sun_Parent, ltime  ; HACK----Not sure if light time is importnant here!!!
      sun2particle = Sun_Parent[0:2,*] + [x, y, z]  ;Calculate the vectors from the sun to the particles
      CSPICE_SPKPOS, Parent_ID, ephemeris_time - REFORM(time), 'J2000', 'NONE', Body, Body2Parent, ltime 
      x_to_Parent = Body2Parent[0,*] - x 
      y_to_Parent = Body2Parent[1,*] - y
      z_to_Parent = Body2Parent[2,*] - z
      dist2Parent = sqrt( (x_to_Parent^2.) + (y_to_Parent^2.) +(z_to_Parent^2.) ) 
      angular_radius_Parent = asin( Parent_radius / dist2Parent)
    endif

; Theta is the Body-Particle-Sun angle (or Sun-Particle-Body angle) for each packet in radians
  sun2particle = sun_particle ; Hack HAck Hack!!! Untested and this could break something comparing Mercury to the Moon!!! Hack HAck Hack!!! 
  theta = Vector_Sep_angle( [x,y,z], sun2particle )

  ; And the Parent-Particle-Sun angle for each packet  
    if keyword_set(Parent_ID) then begin
      ;Parent_theta = Acos( (x_to_Parent*Sun_Parent[0,*] + y_to_Parent*Sun_Parent[1,*] + z_to_Parent*Sun_Parent[2,*]) $
      ;                      / (sqrt((Sun_Parent[0,*]^2.) + (Sun_Parent[1,*]^2.) + (Sun_Parent[2,*]^2.)) * dist2Parent))
    
;      Parent_theta = Acos( (x_to_Parent*sun2particle[0,*] + y_to_Parent*sun2particle[1,*] + z_to_Parent*sun2particle[2,*]) $
;                            / (sqrt((sun2particle[0,*]^2.) + (sun2particle[1,*]^2.) + (sun2particle[2,*]^2.)) * dist2Parent))
      Parent_theta = Vector_Sep_angle( [x_to_Parent,y_to_Parent,z_to_Parent], -sun2particle )
    endif  

; Define an array to hold the fractional illimunation of the number of points for which the illumination fraction is calculated
  illu            = replicate(1., n_elements(x)) 

; Determine the fractional illumination: BODY OCCULTING SUN
  Partial_occult  = where( (theta lt (angular_radius_body + angular_radius_sun))    and (theta gt abs(angular_radius_body - angular_radius_sun)), N_partial_occults )
  total_occult    = where( (theta le abs(angular_radius_body - angular_radius_sun)) and (angular_radius_body ge angular_radius_sun), N_total_occults )
  body_transit    = where( (theta le abs(angular_radius_body - angular_radius_sun)) and (angular_radius_body lt angular_radius_sun), N_transits )

  if N_total_occults gt 0 then illu[total_occult] = 0.                         ; Fully eclipsed particles
  if N_transits gt 0 then illu[body_transit] = ( (angular_radius_sun[body_transit]^2)-(angular_radius_body[body_transit]^2) ) / $
                                                 (angular_radius_sun[body_transit]^2)
  if N_Partial_occults gt 0 then begin
    ; computing the area of overlap in the "lens" region
    ; see http://mathworld.wolfram.com/CircularSegment.html and http://mathworld.wolfram.com/Circle-CircleIntersection.html?affilliate=1 
    ; for Body larger than the Sun. . . (OR FOR ALL CASES?) 
      R_big  = angular_radius_body[Partial_occult]
      r_small = angular_radius_sun[Partial_occult]
      d      = theta[Partial_occult]    
      a      = (1./d)*sqrt(((2.*d*R_big)^2) - ((d^2)-(R_small^2)+(R_big^2))^2) ; this is the chord length in radians connecting the points where the disks intersect   
      d_1    = .5*sqrt(((2.*R_big)^2)-(a^2))                                   ; angular distance from the body center to the "a" chord connecting the points where the disks overlap
      d_2    = d - d_1                                                         ; angular distance from the sun center to the "a" chord connecting the points where the disks overlap
      Area_1 = ((R_big^2)*acos(d_1/R_big)) - d_1*sqrt((R_big^2)-(d_1^2))       ; area of overlap on the side of the Sun, the part of the bodys disk on the far side of the chord. 
      Area_2 = ((r_small^2)*acos(d_2/r_small)) - d_2*sqrt((r_small^2)-(d_2^2)) ; area of overlap on the side of the body, the part of the Suns disk on the far side of the chord. 
      Area_of_overlap           = Area_1 + Area_2   
      illu[Partial_occult] = ((!dpi*(r_small^2)) - Area_of_overlap)/(!dpi*(r_small^2))  ; Iluminated_fraction
  endif
  
; If "Body" is a moon, then also determine the fractional illumination: PARENT BODY OCCULTING SUN
  If keyword_set(Parent_ID) then begin
    partial_Parent_occult  = where( (Parent_theta lt (angular_radius_Parent + angular_radius_sun)) and (Parent_theta gt abs(angular_radius_Parent - angular_radius_sun)), N_partial_Parent_occults )
    total_Parent_occult = where( (Parent_theta le abs(angular_radius_Parent - angular_radius_sun)) and (angular_radius_Parent ge angular_radius_sun), N_total_Parent_occults )
    Parent_transit      = where( (Parent_theta le abs(angular_radius_Parent - angular_radius_sun)) and (angular_radius_Parent lt angular_radius_sun), N_Parent_transits )

    if N_total_Parent_occults gt 0. then illu[total_Parent_occult] = 0. ;Fully eclipsed particles
    if N_Parent_transits gt 0. then illu[Parent_transit] = ((angular_radius_sun[Parent_transit]^2.)-(angular_radius_Parent[Parent_transit]^2.))/(angular_radius_sun[Parent_transit]^2.)
    if N_partial_Parent_occults gt 0. then begin

    ; computing the area of overlap in the "lens" region
    ; see http://mathworld.wolfram.com/CircularSegment.html and http://mathworld.wolfram.com/Circle-CircleIntersection.html?affilliate=1 
    ; Parent larger than the sun, otherwise would meet the above transit criterion. 
      R_big = angular_radius_Parent[partial_Parent_occult]
      r_small = angular_radius_sun[partial_Parent_occult]
      d = Parent_theta[partial_Parent_occult]    
      a = (1./d) * sqrt(((2.*d*R_big)^2.) - ((d^2.)-(R_small^2.)+(R_big^2.))^2.) ; this is the chord length in radians connecting the points where the disks intersect 
      
      d_1 = .5*sqrt(((2.*R_big)^2.)-(a^2.)) ;angular distance from the Parent center to the "a" cord connecting the points where the disks overlap
      d_2 = d - d_1 ;angular distance from the sun center to the "a" cord connecting the points where the disks overlap
      Area_1 = ((R_big^2.)*acos(d_1/R_big)) - d_1*sqrt((R_big^2.)-(d_1^2.)) ;area of overlap on the side of the Sun, the part of the Parents disk on the far side of the chord. 
      Area_2 = ((r_small^2.)*acos(d_2/r_small)) - d_2*sqrt((r_small^2.)-(d_2^2.)) ;area of overlap on the side of the Parent, the part of the Suns disk on the far side of the chord. 
      
      ;sun_mostly_eclipsed = where(d_2 lt 0.) ;if d_2 is negative then the sun is mostly eclipsed
      ;stop
      ;Area_2(sun_mostly_eclipsed) = (!pi*(r_small^2.)) - Area_2(sun_mostly_eclipsed)
      Area_of_overlap = Area_1 + Area_2
      
      ;d_1_B = ((d^2.) - (r_small^2.) + (R_big^2.))/(2.*d)
      
      ;Parent_mostly_eclipsed = where(d_1 le 0.)
      ;sun_mostly_eclipsed = where(d_2 le 0.)
       
      ;Area_of_overlap = (r_small^2.)*acos(d_1/r_small) + (R_big^2.)*acos(d_2/R_big) - $
      ;                  .5*sqrt(((-theta)+r_small+R_big)*(theta+r_small-R_big)*(theta-r_small+R_big)*(theta+r_small+R_big))
                     
      illu[partial_Parent_occult] = ((!dpi*(r_small^2.)) - Area_of_overlap)/(!dpi*(r_small^2.))  ;Iluminated_fraction 0 to 1
    endif 
  endif  ; fraction illumination 0 to 1 due to solar occultation of the parent body (for cases where the 'body' simulated is a moon)

  ; ERROR HANDLING: Make sure no there are no infinite, NaN or negative pixels
  check_bad_pixels = WHERE( FINITE(Illu, /NAN), count_bad_pixels) 
  if count_bad_pixels gt 0 then stop
return, illu
end
