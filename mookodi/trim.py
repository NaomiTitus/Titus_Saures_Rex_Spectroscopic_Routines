from spectroscopic_routines import *

# Trim = True  
# xmin = 0
# xmax = 1024
# ymin = 700
# ymax = 750
# Trim_file_prefix = None
# Trim_ouput_prefix = 't'


# trim(apply_trim = Trim,
# 	files_prefix = None,
# 	file = 'MKD_20211104.0068.fits',
# 	x1 = xmin,
# 	x2 = xmax,
# 	y1 = ymin,
# 	y2 = ymax,
# 	output_prefix = Trim_ouput_prefix)

w_craw, f_craw = np.loadtxt('craw_hi.dat',usecols=(0,1),unpack=True)
# w_lo,f_lo = np.loadtxt('ref_xe_mk_lo.dat',usecols=(0,1),unpack=True)
w,f = np.loadtxt('wtarc.dat',usecols=(0,1),unpack=True)

plt.plot(w_craw,.54*f_craw,label='High res')
# plt.plot(w_lo,f_lo,label='Low res')
plt.plot(w,.1*f,label='IRAF')
plt.legend(loc=1)
plt.show()
