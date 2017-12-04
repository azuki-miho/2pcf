import matplotlib.pyplot as plt
import string
import numpy
import math

ysize = 80
xsize = 80
tpcffile = open("/home/gongjingyu/gcode/astro/2pcf/outcome/20171125/20171125_tpcfv1")
tpcf = numpy.zeros((ysize, xsize))
for i in range(80):
    line = tpcffile.readline()
    linearr = line.strip().split()
    for j in range(80):
        tpcf[i][j] = (math.log(string.atof(linearr[j])+2))

im = plt.imshow(tpcf)
plt.colorbar()
plt.axis("off")
plt.savefig("/home/gongjingyu/gcode/astro/2pcf/outcome/20171125/20171125_tpcfv1_log.eps")
#plt.savefig("./../outcome/20171107_tpcf.eps")
