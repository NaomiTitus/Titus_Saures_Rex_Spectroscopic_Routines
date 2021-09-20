from spectroscopic_routines import *

run_pipeline = True
apply_flux_cal = True

raw_files_prefix = 'a'

grating = 'gr7'
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

IMAGES =  raw_files_prefix+'*43.fits' # or you can a single image
ARC = -1 	# 1: takes arc after, -1: takes arc before or can take file name eg: 'arc.fits'


#############################
#### Trimming Parameters ####
#############################

Trim = True  
xmin = 7
xmax = 1980
ymin = 36
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

construct_flat = False
flat_keyword = 'FLAT' 
flat_files_prefix = 't'
flat_exp_time = 4.0
masterbias = 'Bias.fits'
flat_mode = 'spline'
flat_response = True
flat_display = False
flat_output = 'Flat.fits'


################################
#### Calibration Parameters ####
################################

apply_call = True
files_prefix = None
cal_file_prefix = 't'
masterbias = 'Bias.fits'
masterflat = 'Flat.fits'
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

y0_trace = 27.8
yf_trace = 33.8
manual_trace = False
manual_poly_order = 3
manual_x = [0,500,1175,1976]
manual_y = [22.1,23.87,26.13,27.42]
trace_prominence = 300
tolerance = 3
trace_width = None
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

apwidth = None
skysep = 5
skywidth = 6
skydeg = 0
coaddN = 1
display_apextract = False

##########################################
##### Wavelength Solution Parameters #####
##########################################

reference_spec = 'ref_cuar_gr7_162.dat'
line_list = 'line_list.dat'
wave_min = 3500 
wave_max = 8500 
prominence = 400 
order = 3
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

if (construct_bias is True) or (construct_flat is True): 
	# IMAGES = glob.glob(raw_files_prefix+'*.fits')
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
			display = flat_display,
			output = flat_output)


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

	for k in IMAGES:
	    gr = image_header(k,grating_keyword)
	    if gr == grating:
	    	appy_calibrations(apply_call = apply_call,
	    		identifying_keyword = identifying_keyword,
	    		science_keyword = science_keyword, 
	    		files_prefix = None,
	    		file = cal_file_prefix+k,
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
	    	exposure_type = image_header(k,identifying_keyword)
	    	if exposure_type == science_keyword:
	    		print(k,exposure_type)
	    		number = int(k.split('.')[0].split(raw_files_prefix)[1])
	    		if type(ARC) is int:
	    			arc = 't'+raw_files_prefix+str(number+ARC)+'.fits'
	    		if type(ARC) is str:
	    			arc = 't'+ARC
	    		spec_file = k.split('.fits')[0]
	    		k = 'cfbt'+k
	    		# print (k)
	    		# input()
	    		try:
	    		    image_raw, sky_subtracted, sky, xs, ys, nx, ny, yvals = twodsky(k,object_keyword=object_keyword)
	    		    #
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
	    		    if manual_trace is True:
	    		    	my, myfwhm, trace_c, my_man, trace_c_man = trace_output
	    		    else:
	    		    	my, myfwhm, trace_c = trace_output
	    			#	   
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
	    		except ValueError:
	    		    print ('Fit unsuccessful for '+k+' in twodsky, cannot complete optimal extraction')
	    		    #
	    		    my, myfwhm, trace_c = ap_trace(
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
	    		    	onedspec_optimal = None,
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
	





