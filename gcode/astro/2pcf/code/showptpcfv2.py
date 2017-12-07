import matplotlib.pyplot as plt
import string
import numpy
import math

ysize = 80
xsize = 80
tpcffile = open("/home/gongjingyu/gcode/astro/2pcf/outcome/20171114/20171114_tpcfv1")
tpcf = numpy.zeros((ysize, xsize))
for i in range(80):
    line = tpcffile.readline()
    linearr = line.strip().split()
    for j in range(80):
        tpcf[i][j] = (string.atof(linearr[j]))
        #tpcf[i][j] = (string.atof(linearr[j]))
log_ptpcf = numpy.zeros((xsize/2))
for i in range(40):
    for j in range(80):
        log_ptpcf[i] += tpcf[j][i+40]
    log_ptpcf[i] = math.log(log_ptpcf[i])
x = numpy.linspace(-0.436295,1,40)
plt.plot(x,log_ptpcf)
plt.xlabel(r'$log_{40}r_p/Mpc$')
plt.ylabel(r'$log\xi (projected)$')
plt.savefig("/home/gongjingyu/gcode/astro/2pcf/outcome/20171125/201711_ptpcfv1_log.eps")
#plt.savefig("./../outcome/20171107_tpcf.eps")
