#!/usr/bin/env python

import matplotlib
matplotlib.use('TkAgg')
from pylab import *
from matplotlib.widgets import Slider, Button, RadioButtons
import modetrap_sub
import numpy as np


# Number of beads

nbeads = 1

narg=len(sys.argv)
if narg >=2 :
    nbeads = int(sys.argv[1])

if nbeads == 1:
    ncases = 5
elif nbeads == 2:
    ncases = 3
else:
    print 'Case for',nbeads,'not supported'
    exit()

# This is a wrapper for modetrap that insures that the bead locatations ('xpert')
# are sorted into increasing order.
def mode_wrap(n1,n2,xpert):
    xp=np.array(xpert)

    iargs = np.argsort(xpert)
    xp=xp[iargs]
    res = modetrap_sub.modetrap(n1,n2,xp,.22,.33)
    return res

# Calculate the forward period difference from a set of consecutive periods.
def dp_calc(periods):
    i0 = np.arange(len(periods)-1)
    per = periods[i0]
    dp = periods[i0+1]-periods[i0]
    #return per,dp
    return i0+1,dp


    params = plists[i]
    if nbeads==1:
        periods = mode_wrap(n1,n2,[params[0]])
    
    return per,dp

# Overtone/harmonic values "n" for the first ("n1") and last ("n2") mode to be computed.
n1=1
n2=21
# Default/initial parameters for the location, amplitude (mass), and width of the beads.
locvec0 = [ 0.16, 0.22, 0.33 ]

# Set the initial set of parameters.
loc = []

for i in np.arange(nbeads):
    loc.append(locvec0[i])

# Set the initial x and y limits on the dP vs. P plot
#x1lim = np.pi*(n1-1)
#x2lim = np.pi*(n2)
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
per0, dp0 = dp_calc(exampledata[:,1])
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
h = ((xright-xleft) - (nbeads-1)*xspace)/real(nbeads)
for i in np.arange(nbeads):
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
    locvec = []
    for i in np.arange(nbeads):
        locvec.append(slvec[i].val)
    periods =  mode_wrap(n1,n2,locvec)
    per,dp = dp_calc(periods)
    l.set_xdata(per)
    l.set_ydata(dp)
    draw()

for i in np.arange(nbeads):
    slvec[i].on_changed(update)

show()

