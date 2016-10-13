#This script is for generating Fig.4 in the paper

#Load the pylab and scipy modules
import pylab as pl
import scipy.interpolate as intp

#Generate simulation space and assign material parameters
x = pl.linspace(0e-9, 4000e-9, 1000)
y = pl.linspace(0e-9, 1000e-9, 250)
Ms = 8e5

#Load magnetization and magnetic intensity at t=300 ns and 0 ns
mc = pl.loadtxt('Fig4/alpha-a/m000003.ovf')*Ms
mth = pl.loadtxt('Fig4/alpha-b/m000003.ovf')*Ms
mp1 = pl.loadtxt('Fig4/alpha-c-1/m000003.ovf')*Ms
mp2 = pl.loadtxt('Fig4/alpha-c-2/m000003.ovf')*Ms

mc0 = pl.loadtxt('Fig4/alpha-a/m000000.ovf')*Ms
mth0 = pl.loadtxt('Fig4/alpha-b/m000000.ovf')*Ms
mp10 = pl.loadtxt('Fig4/alpha-c-1/m000000.ovf')*Ms
mp20 = pl.loadtxt('Fig4/alpha-c-2/m000000.ovf')*Ms

bc = pl.loadtxt('Fig4/alpha-a/B_eff000003.ovf')
bth = pl.loadtxt('Fig4/alpha-b/B_eff000003.ovf')
bp1 = pl.loadtxt('Fig4/alpha-c-1/B_eff000003.ovf')
bp2 = pl.loadtxt('Fig4/alpha-c-2/B_eff000003.ovf')

bc0 = pl.loadtxt('Fig4/alpha-a/B_eff000000.ovf')
bth0 = pl.loadtxt('Fig4/alpha-b/B_eff000000.ovf')
bp10 = pl.loadtxt('Fig4/alpha-c-1/B_eff000000.ovf')
bp20 = pl.loadtxt('Fig4/alpha-c-2/B_eff000000.ovf')

print 'Done loading files!'

#Reshape arrays according to simulation space
mc.shape = (y.shape[0], x.shape[0], 3)
mth.shape = (y.shape[0], x.shape[0], 3)
mp1.shape = (y.shape[0], x.shape[0], 3)
mp2.shape = (y.shape[0], x.shape[0], 3)

mc0.shape = (y.shape[0], x.shape[0], 3)
mth0.shape = (y.shape[0], x.shape[0], 3)
mp10.shape = (y.shape[0], x.shape[0], 3)
mp20.shape = (y.shape[0], x.shape[0], 3)

bc.shape = (y.shape[0], x.shape[0], 3)
bth.shape = (y.shape[0], x.shape[0], 3)
bp1.shape = (y.shape[0], x.shape[0], 3)
bp2.shape = (y.shape[0], x.shape[0], 3)

bc0.shape = (y.shape[0], x.shape[0], 3)
bth0.shape = (y.shape[0], x.shape[0], 3)
bp10.shape = (y.shape[0], x.shape[0], 3)
bp20.shape = (y.shape[0], x.shape[0], 3)

#Calculate Energy density for each profile
HBpc = -0.5*pl.add.reduce(mc*bc, axis=2)
HBpth = -0.5*pl.add.reduce(mth*bth, axis=2)
HBpp1 = -0.5*pl.add.reduce(mp1*bp1, axis=2)
HBpp2 = -0.5*pl.add.reduce(mp2*bp2, axis=2)

HBpc0 = -0.5*pl.add.reduce(mc0*bc0, axis=2)
HBpth0 = -0.5*pl.add.reduce(mth0*bth0, axis=2)
HBpp10 = -0.5*pl.add.reduce(mp10*bp10, axis=2)
HBpp20 = -0.5*pl.add.reduce(mp20*bp20, axis=2)

HBpcyavg = pl.average(HBpc, axis=0)
HBpthyavg = pl.average(HBpth, axis=0)
HBpp1yavg = pl.average(HBpp1, axis=0)
HBpp2yavg = pl.average(HBpp2, axis=0)

HBpc0yavg = pl.average(HBpc0, axis=0)
HBpth0yavg = pl.average(HBpth0, axis=0)
HBpp10yavg = pl.average(HBpp10, axis=0)
HBpp20yavg = pl.average(HBpp20, axis=0)

#Remove ground state energy density
HBpcyavg = HBpcyavg - HBpc0yavg
HBpthyavg = HBpthyavg - HBpth0yavg
HBpp1yavg = HBpp1yavg - HBpp10yavg
HBpp2yavg = HBpp2yavg - HBpp20yavg

#Use spline interpolation for getting more intermediate points
# Remoce the first 20 points as they are spurious zeros
xind = pl.where(x*1e9>3800)
xnew = pl.linspace(3.8e-6, 4e-6, 1000)
HBpcnew = intp.spline(x[xind], HBpcyavg[xind], xnew)[19:]
HBpthnew = intp.spline(x[xind], HBpthyavg[xind], xnew)[19:]
HBpp1new = intp.spline(x[xind], HBpp1yavg[xind], xnew)[19:]
HBpp2new = intp.spline(x[xind], HBpp2yavg[xind], xnew)[19:]

#Plot the energy density for the different profiles
fig, ax = pl.subplots()
pl.plot(xnew[19:][::10]*1e6, 10*(pl.log10(HBpcnew)-pl.log10(HBpcnew.max()))[::10], '-', lw=2)
pl.plot(xnew[19:][::10]*1e6, 10*(pl.log10(HBpthnew)-pl.log10(HBpthnew.max()))[::10], 'o', lw=2)
pl.plot(xnew[19:][::10]*1e6, 10*(pl.log10(HBpp1new)-pl.log10(HBpp1new.max()))[::10], '^', lw=2)
pl.plot(xnew[19:][::10]*1e6, 10*(pl.log10(HBpp2new)-pl.log10(HBpp2new.max()))[::10], 'v', lw=2)
pl.xlim(3.8, 3.996)
pl.xlabel(r'$x\,\left(\rm{mm} \right)$', fontsize=24)
pl.ylabel(r'${\cal{E}}\,\left(x,\,\left<y\right>,\,t=300\,\rm{ns} \right)\,\rm{dB}$', fontsize=24)
ax.tick_params(axis='both', which='major', labelsize=18)
ax.tick_params(axis='both', which='minor', labelsize=18)
pl.legend((r'$\alpha_{\rm a}$', r'$\alpha_{\rm b}$', r'$\alpha_{\rm c,\,1}$', r'$\alpha_{\rm c,\,2}$'), loc='best', fontsize=24, frameon=0)
fig.tight_layout()
pl.show()
