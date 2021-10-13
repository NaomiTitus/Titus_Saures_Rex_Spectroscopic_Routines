from spectroscopic_routines import *

run_pipeline = True
apply_flux_cal = False

raw_files_prefix = 'mbxgpP'

grating = 'PG0300'
grating_keyword = 'GRATING'

identifying_keyword = 'OBSTYPE' # identify science / flat / etc
exp_time_keyword = 'EXPTIME'
science_keyword = 'OBJECT'
object_keyword = 'OBJECT'

gain_keyword = 1.507 #'GAINVAL'
readnoise_keyword = 2.225 #'NOISEADU'

spatial_axis = 1
flip_wave_axis = False

#############################
#############################

IMAGES =  raw_files_prefix+'*.fits' # or you can a single image
ARC = 1 	# 1: takes arc after, -1: takes arc before or can take file name eg: 'arc.fits'
final_prefix = 'ct' # this may change depending on the calibration performed


#############################
#### Trimming Parameters ####
#############################

Trim = True  
xmin = 870
xmax = 3030
ymin = 330
ymax = 630
Trim_file_prefix = None
Trim_ouput_prefix = 't'

################################
#### Calibration Parameters ####
################################

apply_call = True,
files_prefix = None
cal_file_prefix = 't'
masterbias = None
masterflat = None
cosmicray_removal = True
contrast = 4
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

y0_trace = 119
yf_trace = 127
trace_prominence = 200
manual_trace = False
manual_poly_order = 3
manual_x = [0,1000,1976]
manual_y = [26.7,29.2,32]
tolerance = 1
trace_width = 'fwhm' # or for a factor of 2 of fwhm: '2-fwhm' or just integers for pixels
poly_order = 3
fmask = (1,) 
nsteps = 20
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

skysep = 5
skywidth = 10
skydeg = 0
display_apextract = False

##########################################
##### Wavelength Solution Parameters #####
##########################################

reference_spec = 'ref_ar_pg300_gr5.dat'
line_list = 'line_list.dat'
wave_min = 3800
wave_max = 9200 
prominence = 10 
order = 2
parameter_file = 'parameters_new'
view_arc = False
display_wave = False

##########################################
######### Flux Calibration ###############
##########################################

display_flux_cal = False

##########################################
############## Pipeline ##################
##########################################


if run_pipeline is True:

	if type(IMAGES) is str:
		IMAGES = glob.glob(IMAGES)

	trim_images = []
	
	for k in IMAGES:
	    gr = image_header(k,grating_keyword)
	    if gr == grating:
	    	trim_images.append(k)
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
	    				print ('ARC length list must match IMAGES length list')
	    				break

	    		spec_file = IMAGES[kk].split('.fits')[0]
	    		k = final_prefix+IMAGES[kk]
	    		# print (k)
	    		# input()
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
	    		my, myfwhm, trace_c = trace_output
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
	    		onedspec, fluxerr, variancespec = ap_extract(
	    			k, 
	    			trace = my, 
	    			apwidth = round(myfwhm), 
	    			skysep = skysep, 
	    			skywidth = skywidth, 
	    			skydeg = skydeg,
	    		    coaddN = coaddN,
	    		    gain_keyword = gain_keyword,
	    		    readnoise_keyword = readnoise_keyword,
	    		    object_keyword = object_keyword,
	    		    display = display_apextract)
	    		#
	    		interact = wavelength(
	    			onedspec,
	    			onedspec_optimal = optimal_spec,
	    			spec_file_name = spec_file,
	    			arc_file = arc,
	    			reference_spec = reference_spec,
	    		    line_list = line_list,
	    		    trace = my,
	    		    trace_fwhm = myfwhm,
	    		    wave_min = wave_min, 
	    		    wave_max = wave_max, 
	    		    prominence = prominence, 
	    		    order = order,
	    		    parameter_file = parameter_file,
	    		    object_keyword = object_keyword,
	    		    flip_wave_axis = flip_wave_axis,
	    		    view_arc = view_arc,
	    		    display = display_wave)
	    		view_arc = interact
	    		    #



if apply_flux_cal is True:
	reduced_data_files = glob.glob('*reduced.dat')
	standard = []
	standard_name = []
	for k in reduced_data_files:
		if k.split('_')[0].lower() == 'ltt3218':
			standard.append(k)
			standard_name.append('ltt3218')
		if k.split('_')[0].lower() == 'eg21':
			standard.append(k)
			standard_name.append('eg21')
		if k.split('_')[0].lower() == 'ltt377':
			standard.append(k)
			standard_name.append('ltt377')
		if k.split('_')[0].lower() == 'feige110':
			standard.append(k)
			standard_name.append('feige110')
		if (k.split('_')[0].lower() == 'cd-32-9927') or (k.split('_')[0].lower() == 'cd-32d9927'):
			standard.append(k)
			standard_name.append('cd-32-9927')
		if k.split('_')[0].lower() == 'ltt7379':
			standard.append(k)
			standard_name.append('ltt7379')
		if k.split('_')[0].lower() == 'ltt7987':
			standard.append(k)
			standard_name.append('ltt7987')
	
	if len(standard) != None:
		for k in reduced_data_files:
			print (k)
			flux_callibration(
		    	standard_reduced_spec = standard[0], 
		    	standard_name = standard_name[0], 
		    	science_spec = k,
		    	display = display_flux_cal)
	else: 
		print('Cannot perform flux callibration')
