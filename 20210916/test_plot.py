import matplotlib.pyplot as plt
import numpy as np

# w, f = np.loadtxt('AT2021uxb_best.dat',skiprows=1,usecols=(0,2),unpack=True)
# w2, f2 = np.loadtxt('AT2021uxb_a2991019_reduced.dat',skiprows=1,usecols=(0,2),unpack=True)

w2, f2 = np.loadtxt('AT2021uxb_a2991019_fluxcal.dat',skiprows=1,usecols=(0,1),unpack=True)

# plt.plot(w,f,label='21')
plt.plot(w2,f2,label='AT2021uxb_a2991019')
plt.legend(loc='best')
plt.xlabel('Wavelength ($\AA$)')
plt.ylabel('$F_\\lambda \\left( sfdsf \\right)$')
plt.show()