#!/bin/py

import base64   #for image to string conversion

def base64png(img) :
  
  with open(img, "rb") as imageFile:
    strimg = base64.b64encode(imageFile.read())
  
  return(strimg)
        
# End
        
