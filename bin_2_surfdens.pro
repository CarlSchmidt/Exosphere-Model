;+
; NAME:
;   Bin_2_SurfDens
;
; PURPOSE:
;   Return the density function (histogram) of latitude and longitude points, with weighting, over a surface grid
;
; CALLING SEQUENCE:
;   Result = Bin_2_SurfDens(V1, V2, W, Radius)
; INPUTS:
;   V1 and V2 = arrays containing the variables.  May be any non-complex
;       numeric type.
;   W = weights to be applied to each element
;
; Keyword Inputs:
;       MIN1:   MIN1 is the minimum V1 value to consider. If this
;               keyword is not specified, then if the smallest value of
;               V1 is greater than zero, then MIN1=0 is used, otherwise
;               the smallest value of V1 is used.
;
;       MIN2:   MIN2 is the minimum V2 value to consider. If this
;               keyword is not specified, then if the smallest value of
;               V2 is greater than zero, then MIN2=0 is used, otherwise
;               the smallest value of V2 is used.
;
;       MAX1:   MAX1 is the maximum V1 value to consider. If this
;               keyword is not specified, then V1 is searched for
;               its largest value.
;
;       MAX2    MAX2 is the maximum V2 value to consider. If this
;               keyword is not specified, then V2 is searched for
;               its largest value.
;
;       BIN1    The size of each bin in the V1 direction (column
;               width).  If this keyword is not specified, the
;               size is set to 1.
;
;       BIN2    The size of each bin in the V2 direction (row
;               height).  If this keyword is not specified, the
;               size is set to 1.
;
; OUTPUTS:
;   The two dimensional density function of the two variables,
;   a longword array of dimensions (m1, m2), where:
;       m1 = Floor((max1-min1)/lon_bin) + 1
;      and  m2 = Floor((max2-min2)/lat_bin) + 1
;   and Result(i,j) is equal to the number of sumultaneous occurences
;   of an element of V1 falling in the ith bin, with the same element
;   of V2 falling in the jth bin. Each occurences is multiplied by its
;   respective weight value W.
;
; RESTRICTIONS:
;   Not usable with complex or string data.
;
; PROCEDURE:
;   Creates a combines array from the two variables, equal to the
;   linear subscript in the resulting 2D histogram, then applies
;   the standard histogram function.
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;   Written by:
;   DMS, Sept, 1992     Written
;   DMS, Oct, 1995      Added MIN, MAX, BIN keywords following
;               suggestion of Kevin Trupie, GSC, NASA/GSFC.
;   CT, RSI, May 2001: Corrected MIN, MAX keywords so that the out-of-range
;               values are ignored rather than truncated to be within range.
;               Allow input arrays with negative values.
; BC, april 2005: added weight
; 
;-
function Bin_2_SurfDens, im1, im2, ww, radius, $
  MIN1 = min1in, MIN2 = min2in, $
  MAX1 = max1in, MAX2 = max2in, $
  BIN1 = lon_bin, BIN2 = lat_bin

  ON_ERROR, 2

  ;Find extents of arrays.
  im1max = MAX(im1, MIN=im1min)
  im2max = MAX(im2, MIN=im2min)

  ;Supply default values for keywords.
  min1 = (N_ELEMENTS(min1in) gt 0) ? min1in : (0 < im1min)
  max1 = (N_ELEMENTS(max1in) gt 0) ? max1in : im1max
  min2 = (N_ELEMENTS(min2in) gt 0) ? min2in : (0 < im2min)
  max2 = (N_ELEMENTS(max2in) gt 0) ? max2in : im2max
  b1 = (N_ELEMENTS(lon_bin) gt 0) ? lon_bin : 1L
  b2 = (N_ELEMENTS(lat_bin) gt 0) ? lat_bin : 1L

  ;Get # of bins for each
  im1bins = FLOOR((max1-min1) / b1) + 1L
  im2bins = FLOOR((max2-min2) / b2) + 1L

  if (im1bins le 0) then MESSAGE, 'Illegal bin size for V1.'
  if (im2bins le 0) then MESSAGE, 'Illegal bin size for V2.'

  noMinTruncation = (min1 eq im1min) and (min2 eq im2min)
  noMaxTruncation = (im1max le max1) and (im2max le max2)
  binSizeOne = (b1 eq 1) and (b2 eq 1)

  im1tmp = im1
  im2tmp = im2
  wwtmp  = ww

  ; Only do the data conversions that are necessary.
  if (min1 ne 0) then im1tmp = TEMPORARY(im1tmp) - min1
  if (min2 ne 0) then im2tmp = TEMPORARY(im2tmp) - min2
  if (b1 ne 1) then im1tmp = TEMPORARY(im1tmp)/b1
  if (b2 ne 1) then im2tmp = TEMPORARY(im2tmp)/b2
  h = im1bins*LONG(TEMPORARY(im2tmp)) + LONG(TEMPORARY(im1tmp))

  ; Construct an array of out-of-range (0) and in-range (1) values.
  in_range = 1
  if (noMinTruncation eq 0) then $ ; set lt min to zero
    in_range = (im1 ge min1) and (im2 ge min2)
  if (noMaxTruncation eq 0) then $ ; set gt max to zero
    in_range = TEMPORARY(in_range) and (im1 le max1) and (im2 le max2)
  ; Set values that are out of range to -1
  h = (TEMPORARY(h) + 1L)*TEMPORARY(in_range) - 1L

  hu = h[uniq(h,sort(h))]
  nhsu  = n_elements(hu)
  hh = fltarr(im1bins, im2bins)

  for ii=1l,nhsu-1l do begin
    ih = hu[ii] mod im1bins
    jh = hu[ii] / im1bins
    hh[ih,jh] = total( ww( where( h eq hu[ii] ) ) )
    ;print,ii,ih,jh,hh[ih,jh]
  endfor
  
  ; Determine the surface area per bin
  quadrangle_area = fltarr(180. / lat_bin)
  for i = -90., 90.-lat_bin, lat_bin do begin
    lon0 = 0.
    lon1 = lon_bin
    lat0 = i
    lat1 = i + lat_bin
    ;quadrangle_area[i] = radius^2 * abs( sin(lat0/!radeg) - sin(lat1/!radeg) ) * abs(lon0/!radeg - lon1/!radeg) / lat_bin
    quadrangle_area[i] = radius^2 * abs( cos(lat0/!radeg) - cos(lat1/!radeg) ) * abs(lon0/!radeg - lon1/!radeg) / lat_bin
  endfor

  surface_area_per_bin = rebin( transpose(quadrangle_area), 360. / lon_bin, 180. / lat_bin)
  hh = hh / surface_area_per_bin
  
  return, hh
end