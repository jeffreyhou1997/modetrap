#!/usr/bin/env python

import matplotlib
matplotlib.use('TkAgg')
from pylab import *
from matplotlib.widgets import Slider, Button, RadioButtons
import modetrap_sub
import numpy as np

#OUTPUT FILE
outdata = []
def output(var1):
    outdata.extend(var1)
    np.savetxt("output.txt",np.array(outdata))

# This is a wrapper for modetrap that insures that the bead locatations ('xpert')
# are sorted into increasing order.
def mode_wrap(n1,n2,xpert):
    xp=np.array(xpert)
    
    iargs = np.argsort(xpert)
    xp=xp[iargs]
    res = modetrap_sub.modetrap(n1,n2,xp,.63,.03)
    return res

# Calculate the forward period difference from a set of consecutive periods.
def dp_calc(periods):
    i0 = np.arange(len(periods)-1)
    per = periods[i0]
    dp = periods[i0+1]-periods[i0]
    #return per,dp
    return i0+1,dp

# Overtone/harmonic values "n" for the first ("n1") and last ("n2") mode to be computed.
n1=1
n2=21
# Default/initial parameters for the location, amplitude (mass), and width of the beads.
locvec0 = [0.0]

# Set the initial set of parameters.
loc = []

loc.append(locvec0[0])

x1lim = n1-0.5
x2lim = n2-0.5
y1lim = 2.0
y2lim = 4.0

xleft = 0.25

fig = figure(1,figsize=(12,6))
ax = subplot(111)
subplots_adjust(left=0.07, bottom=0.22, right=0.975, top=0.95)
periods =  mode_wrap(n1,n2,loc)
per,dp = dp_calc(periods)

l, = plot(per,dp, 'ro-', lw=2, color='red')
exampledata = np.loadtxt('exampledata.txt')
outputdata = np.loadtxt('output.txt')
per0 = exampledata[:,0]
dp0 = exampledata[:,1]
l0, = plot(per0,dp0, 'ko--', lw=3, color='black', markersize=10, mfc='none')

axis([x1lim, x2lim, y1lim, y2lim])
#xlabel(r'$\omega_n\, \rm (frequency)$', fontsize=22)
xlabel(r'$n \,\, \rm (overtone\, number)$', fontsize=22)
ylabel(r'$\omega_{n+1} \,-\, \omega_n$',fontsize=22)

axcolor = 'lightgoldenrodyellow'

axvec = []
slvec = []

xright = 0.93
xleft = 0.1
xspace = 0.1
h = ((xright-xleft) - (1-1)*xspace)/real(1)
for i in np.arange(1):
    xa = xleft + real(i)*(h+xspace)
    axloc  = axes([xa, 0.04, h, 0.06], axisbg=axcolor)
    axs = [ axloc ]
    axvec.append(axs)
    if i==0:
        sloc = Slider(axloc, 'Location ',0.0, 1.0, valinit=locvec0[i])
    
    else:
        sloc = Slider(axloc, '', 0.0, 1.0, valinit=locvec0[i])
    
    slvec.append(sloc)

def update(val):
    var = float(np.round(val*100)/100)
    #print var
    #print outputdata
    l.set_xdata(per)
    var2 = int(np.round(val*100))
    #print var2
    temparray = outputdata[var2*20:var2*20+20]
    l.set_ydata(temparray)
    draw()


slvec[0].on_changed(update)

show()
