#Script to generate Fig. 4 in the paper


#Load the pylab and scipy modules
import pylab as pl
from scipy.optimize import curve_fit

#Generate simulation space
x = pl.linspace(0e-9, 4e-6, 1000)
y = pl.linspace(0e-9, 1000e-9, 250)
xind = pl.intersect1d(pl.where(x*1e9>200)[0], pl.where(x*1e9<3800)[0])
yind = pl.intersect1d(pl.where(y*1e9>200)[0], pl.where(y*1e9<800)[0])

#Load magnetization at t=300 ns
dc = pl.loadtxt('data/normal-incidence/alpha-a/m_y000300.ovf')
dc.shape = (y.shape[0], x.shape[0])
dth = pl.loadtxt('data/normal-incidence/alpha-b/m_y000300.ovf')
dth.shape = (y.shape[0], x.shape[0])
dp1 = pl.loadtxt('data/normal-incidence/alpha-c-1/m_y000300.ovf')
dp1.shape = (y.shape[0], x.shape[0])
dp2 = pl.loadtxt('data/normal-incidence/alpha-c-2/m_y000300.ovf')
dp2.shape = (y.shape[0], x.shape[0])

#Average the magnetization over the width
dcyavg = pl.average(dc[yind], axis=0)
dthyavg = pl.average(dth[yind], axis=0)
dp1yavg = pl.average(dp1[yind], axis=0)
dp2yavg = pl.average(dp2[yind], axis=0)

#Find the zero crossings of the magnetization to find the wavelength and the wavenumber
zero_crossings = pl.where(pl.diff(pl.sign(dcyavg[xind])))[0]
lmbd = 2*pl.average(pl.gradient(x[xind][zero_crossings[2:-3]]))
beta = 2*pl.pi/lmbd

#Find the phase at the ABL interface
#Fit the parameters to the transmission model (refer to the paper)
zcdc = pl.where(pl.diff(pl.sign(dcyavg[xind])))[0]
ldc = x[xind][zcdc][-1]
x0 = 3.8e-6
phdc = (2*pl.pi/lmbd)*(x0-ldc) + pl.pi

def ysig(x, a, G, sigma):
        return a*(pl.exp(-sigma*x))*(pl.cos(beta*x) + G*pl.exp(2*sigma*x)*pl.cos(beta*x+phdc))

poptc, pcovc = curve_fit(ysig, x[xind], dcyavg[xind])

#Repeat the above for the tan hyperbolic profile
zcdth = pl.where(pl.diff(pl.sign(dthyavg[xind])))[0]
ldth = x[xind][zcdth][-1]
phdth = (2*pl.pi/lmbd)*(x0-ldth)

def ysig(x, a, G, sigma):
        return a*(pl.exp(-sigma*x))*(pl.cos(beta*x) + G*pl.exp(2*sigma*x)*pl.cos(beta*x+phdth))

poptth, pcovth = curve_fit(ysig, x[xind], dthyavg[xind])

#Repeat the above for the linear profile
zcdp1 = pl.where(pl.diff(pl.sign(dp1yavg[xind])))[0]
ldp1 = x[xind][zcdp1][-1]
phdp1 = (2*pl.pi/lmbd)*(x0-ldp1)

def ysig(x, a, G, sigma):
        return a*(pl.exp(-sigma*x))*(pl.cos(beta*x) + G*pl.exp(2*sigma*x)*pl.cos(beta*x+phdp1))

poptp1, pcovp1 = curve_fit(ysig, x[xind], dp1yavg[xind])

#Repeat the above for the parabolic profile
zcdp2 = pl.where(pl.diff(pl.sign(dp2yavg[xind])))[0]
ldp2 = x[xind][zcdp2][-1]
phdp2 = (2*pl.pi/lmbd)*(x0-ldp2)

def ysig(x, a, G, sigma):
        return a*(pl.exp(-sigma*x))*(pl.cos(beta*x) + G*pl.exp(2*sigma*x)*pl.cos(beta*x+phdp2))

poptp2, pcovp2 = curve_fit(ysig, x[xind], dp2yavg[xind])

#Print out the loss introduced in the line for the differnt profiles
def alphaf(popt):
        return abs(popt[2])

alphaabs = pl.array([alphaf(poptc), alphaf(poptth), alphaf(poptp1), alphaf(poptp2)])/1e6
print alphaabs

#Print out the return loss for the different profiles
def Gabsf(popt):
        return abs(popt[1])

Gabs = pl.array([Gabsf(poptc), Gabsf(poptth), Gabsf(poptp1), Gabsf(poptp2)])
RL = -20*pl.log10(Gabs)
print RL

def yfit(x, popt):
        return ysig(x, popt[0], popt[1], popt[2])


fig, ax = pl.subplots()
pl.plot(x[xind]*1e6, dp2yavg[xind], 'o-', lw=4)
pl.plot(x[xind]*1e6, yfit(x[xind], poptp2), '-', lw=4)
pl.xlabel(r'$x\,\left(\rm{\mu m} \right)$', fontsize=24)
pl.ylabel(r'${m_{y}}\,\left(x,\,\left<y\right>,\,t=300\,\rm{ns} \right)$', fontsize=24)
pl.legend(('Simulation', 'fit'), loc='best', fontsize=18)
pl.xlim(0.2, 3.8)
ax.tick_params(axis='both', which='major', labelsize=18)
ax.tick_params(axis='both', which='minor', labelsize=18)
fig.tight_layout()
pl.show()
