pro loadgood,red,green,blue	;makes "goodcm"-like color table

;bluei=[1,1,1,0,0,0,0,1]
;greeni=[0,.5,1,1,1,.5,0,0]
;redi=[0,0,0,0,1,1,1,1]

bluei=[1,1,0,0,0,0,0,1]
greeni=[.5,1,1,1,.67,.33,0,0]
redi=[0,0,0,1,1,1,1,1]

minbri=.6	;minimum brightness

;*****************
s=n_elements(bluei)
n=256./s

blue=indgen(256)
green=blue
red=blue

blue=(blue mod n)*s*bluei(blue/n)*(1-minbri)+255*bluei(blue/n)*minbri*(bluei(blue/n) gt 0)
green=(green mod n)*s*greeni(green/n)*(1-minbri)+255*greeni(green/n)*minbri*(greeni(green/n) gt 0)
red=(red mod n)*s*redi(red/n)*(1-minbri)+255*redi(red/n)*minbri*(redi(red/n) gt 0)

red(0)=0
red(255)=255
green(0)=0
green(255)=255
blue(0)=0
blue(255)=255

tvlct,red,green,blue

return

end