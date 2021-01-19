pro Display_Surface_Reimpacts, Release_Flux_map, Reimpacting_Flux_Map, body, ephemeris_time, viewpoint

  ;restore, 'C:\Users\schmidtc\Desktop\remove_me.sav'
  directory = 'C:\IDL\Generic Model V2\read_write\'
  loadct, 13
  lat_lon_bin = 4.

  ;----------------------------------------------------------------------------------------------------------------------------------------------------!
  ;--------------------------------------hack! fix the zero bins...this should be fixed in the scatter_2_sphere binning code!--------------------------!
  ;--------------------------------------then the following fix below wouldn't be needed---------------------------------------------------------------!
  ;----------------------------------------------------------------------------------------------------------------------------------------------------!
    bad = where(Release_Flux_map eq 0, /NULL)
    if (bad NE !NULL) then Release_Flux_map[bad] = (Release_Flux_map[bad - 1] + Release_Flux_map[bad + 1]) / 2.
    bad = where(reimpacting_Flux_map eq 0, /NULL)
    if (bad NE !NULL) then reimpacting_Flux_map[bad] = (reimpacting_Flux_map[bad - 1] + reimpacting_Flux_map[bad + 1]) / 2.
  ;----------------------------------------------------------------------------------------------------------------------------------------------------!
  ;--------------------------------------hack! fix the zero bins...this should be fixed in the scatter_2_sphere code!----------------------------------!
  ;----------------------------------------------------------------------------------------------------------------------------------------------------! 
    
  cspice_bodvrd, Body, 'RADII', 3, Body_radius                         ; Find the simulated body's radius in Km

  ;----------------------------Generate the Latitude & Longitude grid to overlay----------------------------------
  ; Find the angular size
    CSPICE_SPKEZR, body, ephemeris_time, 'J2000', 'LT', viewpoint, BODY_state, ltime
    R_M = 206264.806 * atan(Body_radius[0]/ norm(BODY_state[0:2]))     ; Radius of the body in arcsec

  ; Find the sub-observer planetographic longitude and latitude
    cspice_subpnt, 'Near point: ellipsoid', Body, ephemeris_time, 'IAU_'+body, 'LT+S', viewpoint, spoint, trgepc, srfvec
    f = ( Body_radius[0]-Body_radius[2] ) / Body_radius[0]             ; flatness parameter, mercury is round so this is zero
    cspice_recpgr, body, spoint, Body_radius[0], f, spglon, spglat, spgalt
    SOP_lat = spglat * cspice_dpr()                                     ; planetographic latitude of the sub-observer point on the surface. (0 at equator, 90 at north pole, -90 at south pole)
    SOP_lon = spglon * cspice_dpr()                                     ; planetographic longitude of the sub-observer point on the surface. (0 to 360 west longitude).
    
  ; Find the sub-solar longitude 
    cspice_subslr, 'Near point: ellipsoid', body, ephemeris_time, 'IAU_'+body, 'none', viewpoint, sub_solar_point_planet_frame, trgepc, srfvec ; not sure why we'd need light time
    cspice_recpgr, body, sub_solar_point_planet_frame, Body_radius[0], f, sub_solar_lon, sub_solar_lat, sub_solar_radius ; Convert to planetographic latitude and longitude, in the 'IAU_' body fixed frame
    SSP_Lon = sub_solar_lon * cspice_dpr()  

    ; T3D rotations for spherical mapping use a right-handed coordinate system. Rotation directions are [ positive down, positive right, positive counter-clockwise]
    ; However, longitude increasing w/ time is a left-handed system,
    ; Therefore, we need to flip the maps to be consistent
      RH_Reimpacting_flux_Map = REVERSE( Reimpacting_flux_Map, 1 )
      RH_Release_flux_Map     = REVERSE( Release_flux_Map, 1 )

    ; rebin in lat/lon and make a mesh
      RH_release_flux_map     = rebin(RH_release_flux_map, 360./lat_lon_bin, 180./lat_lon_bin)     ; re-scaling to lat_lon_bin sized degree bins
      RH_Reimpacting_flux_Map = rebin(RH_Reimpacting_flux_Map, 360./lat_lon_bin, 180./lat_lon_bin) ; re-scaling to lat_lon_bin sized degree bins
      map_size             = size(RH_Reimpacting_flux_Map, /dimensions)
      MESH_OBJ, 4, Vertex_List, Polygon_List, REPLICATE(.5, map_size[0], map_size[1]), /closed
      
    ; rotate about the z axis
      quaternion1     = qtcompose([0,0,1], (SOP_lon - 90.)/!radeg)
      SOP_rotated1    = qtvrot(spoint, quaternion1)                           ; sub-observer point
      SSP_rotated1    = qtvrot(sub_solar_point_planet_frame, quaternion1)     ; sub-solar point
    
    ; rotate about the new x axis
      quaternion2     = qtcompose([1,0,0], (SOP_lat - 90.)/!radeg)             
      SOP_rotated2    = qtvrot(SOP_rotated1, quaternion2)                     ; sub-observer point now aligned with Z
      SSP_rotated2    = qtvrot(SSP_rotated1, quaternion2)                     ; sub-solar point has an arbitrary location in the [x,y] plane, we want to align it with negative x
    
    ; using the new x and y of the sub-solar point, find how much to rotate the image by 
    ; in order to align the body-sun vector to the left (in the Y = 0 aka XZ plane, with some negative x component) 
      arctan, SSP_rotated2[0], SSP_rotated2[1], a, a_deg                      ; IDL equivalent to atan2 https://hesperia.gsfc.nasa.gov/ssw/gen/idl/math/arctan.pro
      Sun_to_the_left = 180. + a_deg

  ; Transform the vertices:
    T3D, /RESET      
    T3D, ROTATE = [0.0, 0.0, (SOP_lon - 90.)]                                  ; rotate sub-observer lonfitude about Z, axis facing you "out of the page"
    T3D, ROTATE = [(SOP_lat - 90.), 0.0, 0.]                                   ; rotate sub-observer latitude about X
    T3D, ROTATE = [0.0, 0.0, -Sun_to_the_left]                                 ; rotate one last time about Z to put the sun in the positive x direction
    T3D, TRANSLATE = [0.5, 0.5, 0.5]
    VERTEX_LIST = VERT_T3D(Vertex_List)

  ; Create the window and view:
    WINDOW, 0, XSIZE=512, YSIZE=512
    CREATE_VIEW, WINX=512, WINY=512

  ; Render the mesh:
    scaling = minmax(RH_Reimpacting_flux_Map)
    Reimpacting_Flux_Sphere = POLYSHADE(Vertex_List, Polygon_List, SHADES = bytscl(RH_Reimpacting_flux_Map, scaling[0], scaling[1]), /T3D)
    Release_Flux_Sphere     = POLYSHADE(Vertex_List, Polygon_List, SHADES = bytscl(RH_Release_flux_Map, scaling[0], scaling[1]), /T3D)      
    ;SET_SHADING, LIGHT=[-3000., 0.0, 0.0], REJECT=0 ; it'd be nice if something like this could handle the day/night shadowing, alas i can't get it to work so use white/black contours
    
  ; hack the brightness scaling  
    ratio        = total(RH_Reimpacting_flux_Map) / total(RH_Release_flux_Map)       ; this should be about 1 in the absence of much escape or ionization
    bytscl_ratio = total( bytscl(RH_Reimpacting_flux_Map, scaling[0], scaling[1]) ) / total( bytscl(RH_Release_flux_Map, scaling[0], scaling[1]) )
    HACK_BRIGHTNESS_SCALING = ratio / bytscl_ratio
    
  ; ---------------------lat-lon contours: setup the x-y grid image----------------------------
        s        = size(Reimpacting_Flux_Sphere, /dim)
        xdim     = s[0]
        ydim     = s[1]
        ctr_xpix = xdim / 2.
        ctr_ypix = ydim / 2.
        pix2km   = 2.*Body_radius[0] / ydim                        ; km/pix
        x        = (dindgen(xdim) - ctr_xpix) * pix2km
        y        = (dindgen(ydim) - ctr_ypix) * pix2km
        xsq      = dblarr(xdim, ydim)                              ; xsq is an image the size of the calib-img where the values are the horizontal distance from the "center" in km
        ysq      = dblarr(xdim, ydim)                              ; ysq is an image the size of the calib-img where the values are the vertical distance from the center.
        for j = 0, xdim-1 do ysq[j,*] = y
        for j = 0, ydim-1 do xsq[*,j] = x
    
      ; Get the planetographic lon/lat at each pixel.
      ; This grid looks best if make the arrays big, rotate things to our "Sunward = -X, Observer-Planet = +Z" geometry with interpolation, then shrink shrink them again.
        ysq = rebin(ysq, 1024, 1024)
        xsq = rebin(xsq, 1024, 1024)
        ob     = deprob(xsq, ysq, Body_radius[0], Body_radius[2], SOP_lat, SOP_lon)
        ob.lon[where(ob.lon eq -666, /Null)] = !values.F_nan
        ob.lat[where(ob.lat eq -666, /Null)] = !values.F_nan
        ob.lon = rot(ob.lon, Sun_to_the_left, /interp)             ; rotate things so that the sun is to the left
        ob.lat = rot(ob.lat, Sun_to_the_left, /interp)
        ob.lon = rebin(ob.lon, 1024, 1024)
        ob.lat = rebin(ob.lat, 1024, 1024)
    
      ; Which locations to mark on the lat & lon grid
        lat_contours = indgen(17)*10 - 80                          ; every 10 deg lat
        lon_contours = indgen(12)*30                               ; every 30 deg lon
    
      ; Which countours are on the dayside, and which are on the nightside  
        dayside_range = (SSP_lon + [-90.,90.]) mod 360.
        case 1 of 
          dayside_range[1] lt 90.: day = where( (ob.lon ge dayside_range[0]) or (ob.lon le dayside_range[1]), complement = night)
          dayside_range[0] lt 0. : day = where( (ob.lon ge dayside_range[0] + 360.) or (ob.lon le dayside_range[1]), complement = night)
          (SSP_lon gt 90.) and (SSP_lon lt 270.): day = where( (ob.lon ge dayside_range[0]) and (ob.lon le dayside_range[1]), complement = night)
        endcase
        
        night_lat = ob.lat
        night_lon = ob.lon
        night_lon[day] = !values.F_nan
        night_lat[day] = !values.F_nan

        day_lat = ob.lat
        day_lon = ob.lon
        day_lon[night] = !values.F_nan
        day_lat[night] = !values.F_nan  
        
  cgPS_open, filename = strcompress(directory + 'Surface_Budget.eps'), /ENCAPSULATED, xsize = 6, ysize = 4
    !P.font=1
    !p.charsize = 1.
    !p.charthick = 1.
    device, SET_FONT = 'Helvetica Bold', /TT_FONT

    Pos = cglayout([2,1], xgap = 0, OXMargin=[2,2], OYMargin=[6,0])
    axis_keywords = {xtickformat:'(A1)', ytickformat:'(A1)', ticklen:0.}
    scale_bytscl = minmax(Reimpacting_Flux_Sphere) 
     
    cgimage, float(Release_Flux_Sphere) / HACK_BRIGHTNESS_SCALING, minvalue = scale_bytscl[0], maxvalue = 130, /keep_aspect, /axes, title = 'Release Flux', AXKEYWORDS = axis_keywords, pos = pos[*,0], charsize = 1.3
    ;cgimage, Release_Flux_Sphere / HACK_BRIGHTNESS_SCALING , /keep_aspect, /axes, title = 'Release Flux', AXKEYWORDS = axis_keywords, pos = pos[*,0], charsize = 1.3
    cgcolorbar, range = scaling, title = '[Atoms cm!U-2!N s!U-1!N]', position = [pos[0,0], .17, pos[2,1], .21], charsize = 1.5
    
    cgcontour, Night_lat, /onimage, levels = lat_contours, LABEL = 0, color = 'black', thick = .5
    cgcontour, Night_lon, /onimage, levels = lon_contours, LABEL = 2, color = 'black', thick = .5
    cgcontour, day_lat, /onimage, levels = lat_contours, LABEL = 0, color = 'snow', thick = .5
    cgcontour, day_lon, /onimage, levels = lon_contours, LABEL = 2, color = 'snow', thick = .5
    
    cgimage, Reimpacting_Flux_Sphere, minvalue = scale_bytscl[0], maxvalue = 110, /keep_aspect, /axes, title = 'Reimpacting Flux', pos = pos[*,1], /noerase, AXKEYWORDS = axis_keywords, charsize = 1.3
    
    cgcontour, Night_lat, /onimage, levels = lat_contours, LABEL = 0, color = 'black', thick = .5
    cgcontour, Night_lon, /onimage, levels = lon_contours, LABEL = 2, color = 'black', thick = .5
    cgcontour, day_lat, /onimage, levels = lat_contours, LABEL = 0, color = 'snow', thick = .5
    cgcontour, day_lon, /onimage, levels = lon_contours, LABEL = 2, color = 'snow', thick = .5
  
  cgPS_close

  print, 'sub-observer point planetograpic longitude:', SOP_lon
  print, 'sub-observer point planetograpic latitude: ', SOP_lat
  print, 'sub-solar point longitude:                 ', SSP_lon

  ;  window, 1
  ;      cgimage, Reimpacting_Flux_map, /keep_aspect, /axes
  ;      cgcolorbar, range=float([0,max(Release_flux_Map)]), title = 'Release flux [atoms cm!U-2!N s!U-1!N]', /top, position = [0.125, .89, .925, .92]

  return
end