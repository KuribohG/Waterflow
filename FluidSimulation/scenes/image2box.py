'''scene file type:
# something: comment
box x0 x1 y0 y1 z0 z1: add a box, interval [x0,x1],[y0,y1],[z0,z1]
'''

import cv2
import numpy as np

img = cv2.imread("frog.bmp",0)
h,w = img.shape
print img.shape

fout = open("frog.box","w")

for i in range(h):
	for j in range(w):
		if img[i,j]!=255:
			print >>fout, "box %d %d %d %d %d %d"%(0,10,j,j,w-1-i,w-1-i)
