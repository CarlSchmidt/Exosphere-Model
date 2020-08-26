pro write_avi_animation

;NOTE: To convert postscipt to PNG ImageMagick must be installed on your computer.
;
;See http://www.imagemagick.org/script/index.php




;COMMON Model_shared, Body, Ephemeris_time, Seed, Directory, Particle_data, Line_data, Debug



directory  =         'C:\IDL\Generic Model V2\read_write\' ;Directory where all other files are read from and written to

;Switch into the model's read_write directory   
  CD, Current = Current_Directory   ;Get the current directory
  CD, Directory

  ;Sample an image to get the size:
    cgPS2Raster, strcompress(Directory + 'Output 0_Column_Density.eps'), $
                 strcompress(Directory + 'Output 0_Column_Density.png'), /PNG, density = 600., /portrait
    image = Read_PNG(strcompress(Directory + 'Output 0_Column_Density.png'))
    
    ;check out. . . how's the cropping, resolution ok?
    image_size = size(image,/dimensions)
    window, 0, xs = image_size[1], ys = image_size[2]
    
    ;No? Crop it
    image = image[*,130:1070,0:750]
    cgImage, image 
    image_size = size(image,/dimensions)
    
  ; Set up the video player for output.
     aviFilename = 'Slosh_Column_Density.avi'
     video = IDLffVideoWrite(aviFilename, Format='avi')
     framerate = 10                                                                             ;frames per second
     framedims = [image_size[1], image_size[2]]
     stream = video.AddVideoStream(framedims[0], framedims[1], framerate)
     
     ;add images to the stream
     FOR i = 0, 59 DO BEGIN 
        print, strcompress('Converting image ' + string(i) + ' and addding to AVI stream') 
        cgPS2Raster, strcompress(Directory + 'Output '+ string(i) + '_Column_Density.eps'), $
                     strcompress(Directory + 'Output '+ string(i) + '_Column_Density.png'), /PNG, density = 600., /silent, /portrait
        image = Read_PNG(strcompress(Directory + 'Output '+ string(i) + '_Column_Density.png')) ; Read the high-resolution PNG file.
        image = image[*,130:1070,0:750]
        void = video -> Put(stream, image)                                                      ; Add the high-resolution image to the video stream.
     ENDFOR
    
    ; Clean up by closing the video file and writing its location.
    video -> Cleanup
    Print, 'File "' + aviFilename + '" written to ' + Directory + '.'

;Now do the same for the emission frames

  ; Set up the video player for output.
     aviFilename = 'Slosh_Emission.avi'
     video = IDLffVideoWrite(aviFilename, Format='avi')
     framerate = 10
     framedims = [image_size[1], image_size[2]]
     stream = video.AddVideoStream(framedims[0], framedims[1], framerate)
     
     ;add images to the stream
     FOR i = 0, 59 DO BEGIN 
        print, strcompress('Converting image ' + string(i) + ' and addding to AVI stream') 
        cgPS2Raster, strcompress(Directory + 'Output '+ string(i) + '_Emission.eps'), $
                     strcompress(Directory + 'Output '+ string(i) + '_Emission.png'), /PNG, density = 600., /silent, /portrait
        image = Read_PNG(strcompress(Directory + 'Output '+ string(i) + '_Emission.png')) ; Read the high-resolution PNG file.
        image = image[*,130:1070,0:750]
        void = video -> Put(stream, image)                                                      ; Add the high-resolution image to the video stream.
     ENDFOR
    
    ; Clean up by closing the video file and writing its location.
    video -> Cleanup
    Print, 'File "' + aviFilename + '" written to ' + Directory + '.'
  
  CD, Current_Directory
  
  ;Powerpoint still sometime chokes on AVI movies, have a GIF back up.
  
  
  
  
  
  
end