#Script to generate Fig. 9 in the paper

#Load the pylab and scipy modules
import pylab as pl
import scipy.ndimage

#Create the simulation space
x = pl.linspace(0e-9, 5995e-9, 1200)
y = pl.linspace(0e-9, 2995e-9, 600)

#Load the magnetization at t= 5 ns
mc = pl.loadtxt('data/oblique-incidence/alpha-a/m000500.ovf', usecols=(1,))
mth = pl.loadtxt('data/oblique-incidence/alpha-b/m000500.ovf', usecols=(1,))
mp1 = pl.loadtxt('data/oblique-incidence/alpha-c-1/m000500.ovf', usecols=(1,))
mp2 = pl.loadtxt('data/oblique-incidence/alpha-c-2/m000500.ovf', usecols=(1,))

#Reshape the arrays
mc.shape = (y.shape[0], x.shape[0])
mth.shape = (y.shape[0], x.shape[0])
mp1.shape = (y.shape[0], x.shape[0])
mp2.shape = (y.shape[0], x.shape[0])

num = 1000

#Choose the straight line over which to interpolate
xmin = 792
ymax = 237
xmax = 850
ymin = 66
xline, yline = pl.linspace(xmin, xmax, num), pl.linspace(ymax, ymin, num)
mci = scipy.ndimage.map_coordinates(mc, pl.vstack((yline,xline)))

xmin = 792
ymax = 191
xmax = 835
ymin = 58
xline, yline = pl.linspace(xmin, xmax, num), pl.linspace(ymax, ymin, num)
mthi = scipy.ndimage.map_coordinates(mth, pl.vstack((yline,xline)))

xmin = 794
ymax = 228
xmax = 848
ymin = 67
xline, yline = pl.linspace(xmin, xmax, num), pl.linspace(ymax, ymin, num)
mp1i = scipy.ndimage.map_coordinates(mp1, pl.vstack((yline,xline)))

xmin = 792
ymax = 191
xmax = 835
ymin = 58
xline, yline = pl.linspace(xmin, xmax, num), pl.linspace(ymax, ymin, num)
mp2i = scipy.ndimage.map_coordinates(mp2, pl.vstack((yline,xline)))

#Plot the finlal figure
fig, ax = pl.subplots()
pl.plot(mci[::5], '-', lw=2)
pl.plot(mthi[::5], 'o', lw=2)
pl.plot(mp1i[::5], '^', lw=2)
pl.plot(mp2i[::5], 'v', lw=2)
pl.xlabel('')
pl.ylabel(r'${\cal{E}}\,\left(x,\,\left<y\right>,\,t=300\,\rm{ns} \right)\,\rm{dB}$', fontsize=24)
ax.tick_params(axis='both', which='major', labelsize=18)
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
pl.rc('font', size=18)
xt = [0, 200]
pl.xticks(xt, ['A', 'B'])
pl.legend((r'$\alpha_{\rm a}$', r'$\alpha_{\rm b}$', r'$\alpha_{\rm c,\,1}$', r'$\alpha_{\rm c,\,2}$'), loc='best', fontsize=24, frameon=0)
fig.tight_layout()
pl.show()
