#The Python code to generate Fig.2 in the paper

import pylab as pl

x = pl.linspace(3.8e-6, 3.999e-6, 200)
alphaa = pl.ones(x.shape)*0.1

delalpha = 0.5
sigmax = 40e-9
xprime = 3.9e-6
alphab = delalpha*(1 + pl.tanh((x-xprime)/(sigmax)))

def alphac(x, n):
    a = 1.0/((x.max()-x0)**n)
    return a*((x-x0)**n)

x0 = 3.8e-6

alphac1 = alphac(x, 1.0)
alphac2 = alphac(x, 2.0)

fig, ax = pl.subplots()
pl.plot(x[::2]*1e6, alphaa[::2], '-', lw=2)
pl.plot(x[::2]*1e6, alphab[::2], 'o', lw=2)
pl.plot(x[::2]*1e6, alphac1[::2], '^', lw=2)
pl.plot(x[::2]*1e6, alphac2[::2], 'v', lw=2)
pl.xlim(3.8, 3.996)
pl.xlabel(r'$x\,\left(\rm{\mu m} \right)$', fontsize=24)
pl.ylabel(r'$\alpha\,\left(x \right)$', fontsize=24)
pl.legend((r'$\alpha_{\rm a}$', r'$\alpha_{\rm b}$', r'$\alpha_{\rm c,\,1}$', r'$\alpha_{\rm c,\,2}$'), loc='best', fontsize=24, frameon=0)
ax = pl.gca()
ax.tick_params(axis='both', which='major', labelsize=18)
ax.tick_params(axis='both', which='minor', labelsize=18)
fig.tight_layout()
pl.show()

