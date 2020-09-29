pro scatter_plot_particle_locations,xvals,yvals
  ;;;;;;;;;;
  ;create a plot space into which we will later put observations
  ;;;;;;;;;;

  thetavals= findgen(1000)/1000.*2.0*!pi

  circle_x = cos(thetavals) ;another minus sign because sun is on the left
  circle_y = sin(thetavals)
  plot=plot(circle_x,circle_y,/current,xrange=[3.0,-3],yrange=[-3.,3.], $
    aspect_ratio=1.,xmajor=0,ymajor=0,xminor=0,yminor=0,color=byte(([255,183,76])*0.6), thick=2,$
    xcolor='black',ycolor='black',$
    axis_style=0,antialias=1,/nodata) ;uses /nodata because we actually want circle above observations
    ;we plot again below

  ;plot particle positions
  sym_color='black'
  dot=symbol(xvals,yvals,'o',sym_size=symsize,sym_filled=1,/data,/overplot,sym_color=sym_color)

  ;plot planetary body circle
  circleplot=plot(circle_x,circle_y,/current, $ 
    aspect_ratio=1.,xmajor=0,ymajor=0,xminor=0,yminor=0,color=byte(([255,183,76])),$
    xcolor='white',ycolor='white',thick=2,$
    axis_style=0,antialias=1,/overplot) ;uses /nodata because we actually want circle above observations

  return
end