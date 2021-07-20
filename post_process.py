from spectroscopic_routines import *

# Combining spectra

combine_dat(['MLTJ104704.69-352725.3_a2661022_reduced.dat','MLTJ104704.69-352725.3_a2661024_reduced.dat'],\
	out_file = 'MLTJ104704.69-352725.3_combined',\
	spectrum_keyword='spectrum_optimal',\
	wavelength_keyword='wavelength',\
	combine_method='sum',\
	display=True,\
	sigma=3,accuray=2)


# Normalisation

normalise(file = 'MLTJ104704.69-352725.3_combined.dat',
	wavelength_keyword = 'median waves',
	spectrum_keyword = 'clipped sum fluxes',
	exclude_range = None, #  [4354,4540,5786,5989] 4354 4540 5786 5989
	interactive = False,
	display = True,
	out_file = 'MLTJ104704.69-352725.3',
	poly_order=8,)


# Smoothing

smoothing(file = 'MLTJ104704.69-352725.3_combined.dat',
	wavelength_keyword = 'median waves',
	spectrum_keyword = 'clipped sum fluxes',
	window = 3,
	sigma = None, #10,
	out_file = 'MLTJ104704.69-352725.3',
	display = False)