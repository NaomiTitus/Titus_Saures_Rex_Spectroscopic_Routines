from spectroscopic_routines import *

run_pipeline = True
apply_flux_cal = False

raw_files_prefix = 'a'

grating = 'gr6'
grating_keyword = 'GRATING'

identifying_keyword = 'EXPTYPE' # identify science / flat / etc
exp_time_keyword = 'EXPTIME'
science_keyword = 'SCIENCE'
object_keyword='OBJECT'
gain_keyword = 'GAINVAL'
readnoise_keyword = 'NOISEADU'

spatial_axis = 1
flip_wave_axis = True


#############################
#############################

IMAGES = raw_files_prefix+'*.fits' # or you can a single image or a list of IMAGES
ARC = -1 	# 1: takes arc after, -1: takes arc before or can take file name eg: 'arc.fits'
final_prefix = 'cfbt' # this may change depending on the calibration performed


#############################
#### Trimming Parameters ####
#############################

Trim = True  
xmin = 3
xmax = 1980
ymin = 40
ymax = 98
Trim_file_prefix = None
Trim_ouput_prefix = 't'

#############################
## Master Bias Parameters ###
#############################

construct_bias = False
bias_keyword = 'BIAS'
bias_files_prefix = 't'
bias_output = 'Bias.fits'

#############################
## Master Flat Parameters ###
#############################

construct_flat = True
flat_keyword = 'FLAT' 
flat_files_prefix = 't'
flat_exp_time = 8.0
masterbias = 'Bias.fits'
flat_mode = 'spline' # or poly
flat_response = True
sp_ext = 0				#spline ext
sp_k = 2  				#spline k
sp_s = 0.001 			#spline s
flat_poly_order = 2 	#if poly is chosen, you must provide order
flat_display = False
flat_output = 'Flat.fits'


################################
#### Calibration Parameters ####
################################

apply_call = True
files_prefix = None
cal_file_prefix = 't'
masterbias = 'bBias.fits'
masterflat = 'bFlat.fits'
cosmicray_removal = True
contrast = 6
cr_threshold = 8
neighbor_threshold = 5
error = None
mask = None
background = None
maxiter = 4
border_mode = 'mirror'

################################
###### Trace Parameters ########
################################

y0_trace = 23.7
yf_trace = 28.7
manual_trace = False
manual_poly_order = 3
manual_x = [3,526,1383,1972]
manual_y = [29.75,32,35,36]
trace_prominence = 300
tolerance = 3
trace_width = '3-fwhm' # or for a factor of 2 of fwhm: '2-fwhm' or just integers for pixels
poly_order = 2
fmask = (1,) 
nsteps = 25
recenter = True 
prevtrace = (0,)
bigbox = 15 
Saxis = spatial_axis
display_trace = False

################################
## Optimal Extract Parameters ##
################################

display_optimal = False

################################
##### Apextract Parameters #####
################################

column = None # When None will choose central cross section pixel, otherwise specify pixel number
skysep = 5
skywidth = 6
skydeg = 0
display_apextract = False
interact = False

##########################################
##### Wavelength Solution Parameters #####
##########################################

reference_spec = 'ref_cuar_gr6_118.dat'
line_list = 'line_list.dat'
wave_min = 4000 
wave_max = 7000 
prominence = 100 
order = 4
std_factor = 2 # STD factor to exclude matched lines
parameter_file = 'parameters_new'
view_arc = False
display_wave = True

##########################################
######### Flux Calibration ###############
##########################################

display_flux_cal = False
science_spec = '*_reduced.dat' # Or provide a list: ['sfds.dat','sfsdf.dat'] or single file: asdfds.dat
standard_reduced_spec = 'CD-32d9927_a2781019_reduced.dat'
standard_file = 'fcd32d9927.dat'
standard_name = 'cd-32-9927'

# https://ftp.eso.org/pub/usg/standards/ctiostan/ --> Download fluxcal standards here.

##########################################
############## Pipeline ##################
##########################################

if (construct_bias is True) or (construct_flat is True): 
	IMAGES = glob.glob(raw_files_prefix+'*.fits')
	bias_files = []
	flat_files = []
	for k in IMAGES:
	    gr = image_header(k,grating_keyword)
	    if (gr == grating):
	    	if  (image_header(k,identifying_keyword) == 'BIAS') or (image_header(k,identifying_keyword) == 'FLAT'):
	    		trim(
					apply_trim = Trim,
				    files_prefix = None,
				    file = k,
				    x1 = xmin,
				    x2 = xmax,
				    y1 = ymin,
				    y2 = ymax,
				    output_prefix = Trim_ouput_prefix)
	    		#
	    	if  (image_header(k,identifying_keyword) == 'BIAS'):
	    		bias_files.append(k)
	    	if  (image_header(k,identifying_keyword) == 'FLAT'):
	    		flat_files.append(k)
	if len(bias_files) != 0:
	    master_bias(
			construct_bias = construct_bias,
			identifying_keyword = identifying_keyword,
			bias_keyword = bias_keyword, 
			files_prefix = bias_files_prefix,
			raw_prefix = raw_files_prefix,
			output = bias_output)
	    		#
	if len(flat_files) != 0:
	    master_flat(
			construct_flat = construct_flat,
			identifying_keyword = identifying_keyword,
			flat_keyword = flat_keyword, 
			files_prefix = flat_files_prefix,
			raw_prefix = raw_files_prefix,
			exp_time_keyword = exp_time_keyword,
			exp_time = flat_exp_time, 
			masterbias = masterbias, 
			spatial_axis = spatial_axis,
			mode = flat_mode,
			response = flat_response,
			sp_ext = sp_ext,
			sp_k = sp_k,
			sp_s = sp_s, 
			flat_poly = flat_poly_order, 
			display = flat_display,
			output = flat_output)


if run_pipeline is True:


	if (type(IMAGES) is str):
		Im = glob.glob(IMAGES)

		IMAGES = []
		for k in Im:
			if (image_header(k,identifying_keyword) == science_keyword) and \
			(image_header(k,grating_keyword) == grating):
				IMAGES.append(k)
	elif type(IMAGES) is list:
		IMAGES = IMAGES


	ARC_IMAGES = []

	if type(ARC) is int:
		for i in IMAGES:
			if (image_header(i,identifying_keyword) == science_keyword) and \
			(image_header(i,grating_keyword) == grating):
				number = int(i.split('.')[0].split(raw_files_prefix)[1])
				ARC_IMAGES.append(raw_files_prefix+str(number+ARC)+'.fits')
	if type(ARC) is str:
		ARC_IMAGES.append(ARC)
	if type(ARC) is list:
		for i in ARC:
			ARC_IMAGES.append(i)

	IMAGES = list(set(IMAGES))

	trim_images = []
	
	for k,j in zip(IMAGES,ARC_IMAGES):
	    gr = image_header(k,grating_keyword)
	    if gr == grating:
	    	print (k)
	    	trim_images.append(k)
	    gr = image_header(j,grating_keyword)
	    if gr == grating:
	    	trim_images.append(j)

	#
	for k in trim_images:
	    trim(
			apply_trim = Trim,
	    	files_prefix = None,
	    	file = k,
	    	x1 = xmin,
	    	x2 = xmax,
	    	y1 = ymin,
	    	y2 = ymax,
	    	output_prefix = Trim_ouput_prefix)

	for kk in range(len(IMAGES)):
	    gr = image_header(IMAGES[kk],grating_keyword)
	    if gr == grating:
	    	appy_calibrations(apply_call = apply_call,
	    		identifying_keyword = identifying_keyword,
	    		science_keyword = science_keyword, 
	    		files_prefix = None,
	    		file = cal_file_prefix+IMAGES[kk],
	    		masterbias = masterbias,
	    		masterflat = masterflat,
	    		cosmicray_removal = cosmicray_removal,
	    		contrast = contrast,
	    		cr_threshold = cr_threshold,
	    		neighbor_threshold = neighbor_threshold,
	    		error = error,
	    		mask = mask,
	    		background = background, 
	    		gain_keyword = gain_keyword,
	    		readnoise_keyword = readnoise_keyword, 
	    		maxiter = maxiter, 
	    		border_mode = border_mode)
	    	#
	    	exposure_type = image_header(IMAGES[kk],identifying_keyword)
	    	if exposure_type == science_keyword:
	    		print(IMAGES[kk],exposure_type)
	    		number = int(IMAGES[kk].split('.')[0].split(raw_files_prefix)[1])
	    		if type(ARC) is int:
	    			arc = 't' + raw_files_prefix+str(number+ARC)+'.fits'
	    		if type(ARC) is str:
	    			arc = 't' + ARC
	    		if type(ARC) is list:
	    			if len(ARC) == len(IMAGES):
	    				arc = 't' + ARC[kk]
	    			else:
	    				print ('Arc length list must match Images length list')
	    				break

	    		spec_file = IMAGES[kk].split('.fits')[0]
	    		k = final_prefix+IMAGES[kk]
	    	
	    		try:
	    		    image_raw, sky_subtracted, sky, xs, ys, nx, ny, yvals = twodsky(k,object_keyword=object_keyword)
	    		    optimal_spec = True
	    		    #
	    		except ValueError:
	    		    print ('Fit unsuccessful for '+k+' in twodsky, cannot complete optimal extraction')
	    		    optimal_spec = None

	    		trace_output = ap_trace(
	    			k,
	    			y0_trace = y0_trace,
	    			yf_trace = yf_trace, 
	    			trace_prominence = trace_prominence,
	    			manual_trace = manual_trace,
	    			manual_poly_order = manual_poly_order,
	    			manual_x = manual_x,
	    			manual_y = manual_y,
	    			tolerance = tolerance,
	    			trace_width = trace_width,
	    			poly_order = poly_order, 
	    		    object_keyword = object_keyword, 
	    		    fmask = fmask, 
	    		    nsteps = nsteps,
	    		    recenter = recenter, 
	    		    prevtrace = prevtrace, 
	    		    bigbox = bigbox, 
	    		    Saxis = Saxis, 
	    		    display = display_trace)
	    		my, myfwhm, trace_c, poly = trace_output
	    		if optimal_spec is True:
	    			try:    
	    				optimal_spec = optimal(
	    		    		image_raw, 
	    		    		sky = sky, 
	    		    		xs = xs, 
	    		    		ys = ys, 
	    		    		nx = nx, 
	    		    		ny = ny, 
	    		    		yvals = yvals,
	    		    		trace_c = trace_c,
	    		    		display = display_optimal)
	    			except:
	    				optimal_spec = None

	    		#
	    		onedspec, fluxerr, variancespec, snr_spec = ap_extract(
	    			k, 
	    			trace = my,
	    			poly = poly, 
	    			apwidth = round(myfwhm), 
	    			skysep = skysep, 
	    			skywidth = skywidth, 
	    			skydeg = skydeg,
	    			column = column,
	    		    gain_keyword = gain_keyword,
	    		    readnoise_keyword = readnoise_keyword,
	    		    object_keyword = object_keyword,
	    		    spatial_axis = spatial_axis,
	    		    display = display_apextract,
	    		    interact = interact)
	    		#
	    		interact = wavelength(
	    			onedspec,
	    			onedspec_optimal = optimal_spec,
	    			spec_file_name = spec_file,
	    			snr_spec = snr_spec,
	    			arc_file = arc,
	    			reference_spec = reference_spec,
	    		    line_list = line_list,
	    		    trace = my,
	    		    trace_fwhm = myfwhm,
	    		    wave_min = wave_min, 
	    		    wave_max = wave_max, 
	    		    prominence = prominence, 
	    		    order = order,
	    		    std_factor = std_factor,
	    		    parameter_file = parameter_file,
	    		    object_keyword = object_keyword,
	    		    flip_wave_axis = flip_wave_axis,
	    		    view_arc = view_arc,
	    		    display = display_wave)
	    		view_arc = interact
	    		    #



# if apply_flux_cal is True:
# 	reduced_data_files = glob.glob('*reduced.dat')
# 	standard = []
# 	standard_name = []
# 	for k in reduced_data_files:
# 		if k.split('_')[0].lower() == 'ltt3218':
# 			standard.append(k)
# 			standard_name.append('ltt3218')
# 		if k.split('_')[0].lower() == 'eg21':
# 			standard.append(k)
# 			standard_name.append('eg21')
# 		if k.split('_')[0].lower() == 'ltt377':
# 			standard.append(k)
# 			standard_name.append('ltt377')
# 		if k.split('_')[0].lower() == 'feige110':
# 			standard.append(k)
# 			standard_name.append('feige110')
# 		if (k.split('_')[0].lower() == 'cd-32-9927') or (k.split('_')[0].lower() == 'cd-32d9927'):
# 			standard.append(k)
# 			standard_name.append('cd-32-9927')
# 		if k.split('_')[0].lower() == 'ltt7379':
# 			standard.append(k)
# 			standard_name.append('ltt7379')
# 		if k.split('_')[0].lower() == 'ltt7987':
# 			standard.append(k)
# 			standard_name.append('ltt7987')
	
# 	if len(standard) != None:
# 		for k in reduced_data_files:
# 			print (k)
# 			flux_callibration(
# 		    	standard_reduced_spec = standard[0], 
# 		    	standard_name = standard_name[0], 
# 		    	science_spec = k,
# 		    	display = display_flux_cal)
# 	else: 
# 		print('Cannot perform flux callibration')


# Flux calibration

if apply_flux_cal is True:
	if (type(science_spec) is str) and ('*' in science_spec):
		files =  glob.glob(science_spec)
	elif type(science_spec) is list:
		files = science_spec
	elif (type(science_spec) is str) and ('*' not in science_spec):
		files =  [science_spec]
	
	for k in files:
		flux_callibration(
			science_spec = k,
			standard_reduced_spec = standard_reduced_spec,
			standard_name = standard_name,
			standard_file = standard_file,
			display = display_flux_cal)






