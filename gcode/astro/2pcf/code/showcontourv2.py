import matplotlib.pyplot as plt
import string
import numpy
import math

ysize = 80
xsize = 80
x = numpy.linspace(-40, 40, xsize)
y = numpy.linspace(-40, 40, ysize)
X, Y = numpy.meshgrid(x, y)
tpcffile = open("/home/gongjingyu/gcode/astro/2pcf/outcome/20171125/20171125_tpcfv1tolinear")
tpcf = numpy.zeros((ysize, xsize))
for i in range(80):
    line = tpcffile.readline()
    linearr = line.strip().split()
    for j in range(80):
        tpcf[i][j] = (math.log(string.atof(linearr[j])+2))
        #tpcf[i][j] = (string.atof(linearr[j]))

plt.contourf(X, Y, tpcf, 7, alpha=0)
C = plt.contour(X, Y, tpcf, 7, colors='black', linewidth=0.5)
plt.clabel(C, inline=True, fontsize=10)
plt.savefig("/home/gongjingyu/gcode/astro/2pcf/outcome//20171125/20171125_tpcfv1tolinearlog_contour.eps")
#plt.savefig("./../outcome/20171107_tpcf.eps")
