# python modules that we will use
import os
import numpy as np
import numpy.ma as ma
from astropy.io import fits
import glob

from astropy.modeling import models
import astropy.units as u
from astropy.convolution import convolve, Box1DKernel
from astropy.stats import gaussian_sigma_to_fwhm
from astropy.stats import sigma_clip
from specutils import Spectrum1D, SpectralRegion

import matplotlib
import matplotlib.pylab as plt
import matplotlib.cm as cm
import lacosmic

from specutils.spectra import Spectrum1D
from specutils.fitting import fit_lines
from specutils.manipulation import (box_smooth, gaussian_smooth, trapezoid_smooth)

from scipy.signal import find_peaks, peak_widths, peak_prominences
import scipy.signal
from scipy.optimize import curve_fit, fmin
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import SmoothBivariateSpline
from scipy.interpolate import interp1d, LSQUnivariateSpline, LSQBivariateSpline
from scipy import interpolate
from specutils.analysis import gaussian_sigma_width, gaussian_fwhm
from scipy.stats import linregress

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def show_image(image, lower=-1, upper=3, extent=None):
    sample = sigma_clip(image)
    vmin = sample.mean() + lower * sample.std()
    vmax = sample.mean() + upper * sample.std()
    plt.figure(figsize=(15, 7))
    plt.imshow(image, origin='lower', cmap='gray', aspect='auto', vmin=vmin, vmax=vmax, extent=extent)
    plt.xlabel('Column Number')
    plt.ylabel('Row Number')

def image_header(image,keyword):
    header = fits.getheader(image)
    image = fits.getdata(image)
    #image, header = fits.getdata(image, header=True)
    #print (header)
    param = header[keyword]
    return param



def _gaus(x, a, b, x0, sigma):
    """
    Simple Gaussian function, for internal use only
    Parameters
    ----------
    x : float or 1-d numpy array
        The data to evaluate the Gaussian over
    a : float
        the amplitude
    b : float
        the constant offset
    x0 : float
        the center of the Gaussian
    sigma : float
        the width of the Gaussian
    Returns
    -------
    Array or float of same type as input (x).
    """
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2)) + b

def trim(apply_trim,files_prefix,file,x1,x2,y1,y2,output_prefix):
    if apply_trim is not None or not False:
        cwd = os.getcwd()
        data_dir = cwd
        if files_prefix is not None and file is None:
            files = glob.glob(files_prefix+'*.fits')
            for i in files:
                header = fits.getheader(os.path.join(data_dir,i))
                image = fits.getdata(os.path.join(data_dir,i))
                # image, header = fits.getdata(os.path.join(data_dir,i), header=True)
                trim = image[y1:y2,x1:x2]
                fits.writeto(filename=output_prefix+i,data=np.array(trim),header=header,overwrite=True)
        if file is not None:
            header = fits.getheader(file)
            image = fits.getdata(file)
            # image, header = fits.getdata(file, header=True)
            trim = image[y1:y2,x1:x2]
            fits.writeto(filename=output_prefix+file,data=np.array(trim),header=header,overwrite=True)


# trim(files_prefix='a*',x1=3,x2=1980,y1=40,y2=98,output_prefix='t')

# input()


def master_bias(construct_bias,identifying_keyword,bias_keyword,files_prefix,raw_prefix,output):
    if construct_bias is True:
        k = 0
        files = glob.glob(files_prefix+raw_prefix+'*.fits')
        for i in files:    
            image, header = fits.getdata(i, header=True)
            if header[identifying_keyword] == bias_keyword:
                if k == 0:
                    all_data = image
                    k += 1
                elif k > 0:
                    all_data = np.dstack( (all_data, image) )
                    k += 1
        if k > 1:
            # do median across whole stack
            bias = np.nanmedian(all_data, axis=2)
            hduOut = fits.PrimaryHDU(bias)
            hduOut.writeto(output, overwrite=True)
        else:
            print('No bias files')

# master_bias(identifying_keyword ='EXPTYPE',bias_keyword='BIAS', files_prefix='t',output='Bias.fits')


def master_flat(construct_flat,identifying_keyword,flat_keyword,files_prefix,raw_prefix,exp_time_keyword,exp_time,masterbias,
    output,spatial_axis,mode,response,display):
    if construct_flat is True:
        k = 0
        bias = fits.getdata(masterbias)
        # cwd = os.getcwd()
        # data_dir = cwd
        files = glob.glob(files_prefix+raw_prefix+'*.fits')
        for i in files:    
            image, header = fits.getdata(i, header=True)
            if (header[identifying_keyword] == flat_keyword) and (header[exp_time_keyword] == str(exp_time)):
                print (i)
                b_image = image - bias
                # check for bad regions (not illuminated) in the spatial direction
                ycomp = b_image.sum(axis=spatial_axis) # compress to spatial axis only
                illum_thresh = 0.8 # value compressed data must reach to be used for flat normalization
                ok = np.where( (ycomp>= np.nanmedian(ycomp)*illum_thresh) )
                #
                # assume a median scaling for each flat to account for possible different exposure times
                if (k == 0):
                    all_data = b_image / np.nanmedian(b_image[ok,:])
                    k +=1
                elif (k > 0):
                    all_data = np.dstack( (all_data, b_image / np.nanmedian(b_image[ok,:])) )
                    k += 1
        #
        if k > 1:
            # do median across whole stack of flat images
            flat_stack = np.nanmedian(all_data, axis=2)
            # define the wavelength axis
            Waxis = 0
            # add a switch in case the spatial/wavelength axis is swapped
            if spatial_axis == 0:
                Waxis = 1
            #
            if response is True:
                xdata = np.arange(all_data.shape[1]) # x pixels
                # sum along spatial axis, smooth w/ 5pixel boxcar, take log of summed flux
                print (len(flat_stack[ok].mean(axis=Waxis)))
                flat_1d = np.log10(convolve(flat_stack[ok].mean(axis=Waxis), Box1DKernel(5)))
                #
                if mode=='spline':
                    spl = UnivariateSpline(xdata, flat_1d, ext=0, k=2 ,s=0.001)
                    flat_curve = 10.0**spl(xdata)
                elif mode=='poly':
                    # fit log flux with polynomial
                    flat_fit = np.polyfit(xdata, flat_1d, flat_poly)
                    # get rid of log
                    flat_curve = 10.0**np.polyval(flat_fit, xdata)
                #
                if display is True:
                    plt.figure()
                    plt.plot(10.0**flat_1d)
                    plt.plot(xdata, flat_curve,'r')
                    plt.show()
                #
                # divide median stacked flat by this RESPONSE curve
                flat = np.zeros_like(flat_stack)
                #
                if spatial_axis == 1:
                    for i in range(flat_stack.shape[Waxis]):
                        flat[i,:] = flat_stack[i,:] / flat_curve
                else:
                    for i in range(flat_stack.shape[Waxis]):
                        flat[:,i] = flat_stack[:,i] / flat_curve
            else:
                flat = flat_stack
            #
            # if display is True:
            #       plt.figure()
            #       plt.imshow(flat, origin='lower',aspect='auto')
            #       plt.show()
            #
            # write output to disk for later use
            hduOut = fits.PrimaryHDU(flat)
            hduOut.writeto(output, overwrite=True)
            #
            return flat ,ok[0]
            #
        else:
            print('No flat files')

def appy_calibrations(apply_call,identifying_keyword,science_keyword,files_prefix,file,masterbias,masterflat,cosmicray_removal
    ,contrast,cr_threshold,neighbor_threshold,error,mask,background, gain_keyword, readnoise_keyword,maxiter, border_mode):
    if apply_call is not None or False:
        if (masterbias != None) or (masterflat != None):
            bias = fits.getdata(masterbias)
            flat = fits.getdata(masterflat)
        # cwd = os.getcwd()
        # data_dir = cwd
        if files_prefix is not None:
            files = glob.glob(files_prefix+'*.fits')
            for i in files:  
                header = fits.getheader(i)
                image = fits.getdata(i)  
                # image, header = fits.getdata(i, header=True)
                if header[identifying_keyword] == science_keyword:
                    if type(gain_keyword) == str:
                        effective_gain = float(header[gain_keyword])  #1.145
                    elif type(gain_keyword) == float:
                        effective_gain = gain_keyword
                    if type(readnoise_keyword) == str:
                        readnoise = float(header[readnoise_keyword]) 
                    elif type(readnoise_keyword) == float:
                        readnoise = readnoise_keyword
                    if (header[identifying_keyword] == science_keyword):
                        if (masterbias == None) or (masterflat == None):
                            c_image = lacosmic.lacosmic(image, contrast=contrast, cr_threshold=cr_threshold, neighbor_threshold=neighbor_threshold, 
                                error=error, mask=mask, background=background, effective_gain=effective_gain, readnoise= readnoise,
                                 maxiter=maxiter, border_mode=border_mode)
                            fits.writeto(filename='c'+i,data=np.array(c_image[0]),header=header,overwrite=True)
                        else:
                            fb_image = (image - bias)/flat
                            if cosmicray_removal == False:
                                fits.writeto(filename='fb'+i,data=np.array(fb_image),header=header,overwrite=True)
                            elif cosmicray_removal == True:
                                # See python's lacosmic routine for more options
                                cfb_image = lacosmic.lacosmic(fb_image, contrast=contrast, cr_threshold=cr_threshold, neighbor_threshold=neighbor_threshold, 
                                    error=error, mask=mask, background=background, effective_gain=effective_gain, readnoise= readnoise,
                                     maxiter=maxiter, border_mode=border_mode)
                                # print ('ok')
                                fits.writeto(filename='fb'+i,data=np.array(fb_image),header=header,overwrite=True)
                                fits.writeto(filename='cfb'+i,data=np.array(cfb_image[0]),header=header,overwrite=True)
        if file is not None or not False:
            header = fits.getheader(file)
            image = fits.getdata(file)
            if header[identifying_keyword] == science_keyword:
                if type(gain_keyword) == str:
                        effective_gain = float(header[gain_keyword])  #1.145
                elif type(gain_keyword) == float:
                    effective_gain = gain_keyword
                if type(readnoise_keyword) == str:
                    readnoise = float(header[readnoise_keyword]) 
                elif type(readnoise_keyword) == float:
                    readnoise = readnoise_keyword
                if (header[identifying_keyword] == science_keyword):
                    if (masterbias == None) or (masterflat == None):
                            c_image = lacosmic.lacosmic(image, contrast=contrast, cr_threshold=cr_threshold, neighbor_threshold=neighbor_threshold, 
                                error=error, mask=mask, background=background, effective_gain=effective_gain, readnoise= readnoise,
                                 maxiter=maxiter, border_mode=border_mode)
                            fits.writeto(filename='c'+file,data=np.array(c_image[0]),header=header,overwrite=True)
                    else:
                        fb_image = (image - bias)/flat
                        if cosmicray_removal == False:
                            fits.writeto(filename='fb'+file,data=np.array(fb_image),header=header,overwrite=True)
                        elif cosmicray_removal == True:
                            # See python's lacosmic routine for more options
                            cfb_image = lacosmic.lacosmic(fb_image, contrast=contrast, cr_threshold=cr_threshold, neighbor_threshold=neighbor_threshold, 
                                error=error, mask=mask, background=background, effective_gain=effective_gain, readnoise= readnoise,
                                 maxiter=maxiter, border_mode=border_mode)
                            # print ('ok')
                            fits.writeto(filename='fb'+file,data=np.array(fb_image),header=header,overwrite=True)
                            fits.writeto(filename='cfb'+file,data=np.array(cfb_image[0]),header=header,overwrite=True)



# appy_callibrations(identifying_keyword ='EXPTYPE',science_keyword='SCIENCE', files_prefix='t',
#     masterbias='Bias.fits',masterflat='Flat.fits',cosmicray_removal=True,contrast=6,cr_threshold=8,
#     neighbor_threshold=5,error=None,mask=None,background=None, gain_keyword='GAINVAL',readnoise_keyword='NOISEADU',
#     maxiter=4, border_mode='mirror')

def twodsky(image,object_keyword):
    # cwd = os.getcwd()
    # data_dir = cwd

    image_name = image.split('.')[0]+'.fits'
    image, header = fits.getdata(image, header=True)   
    target_name = header[object_keyword]

    ny, nx = image.shape
    cy, cx = ny//2, nx//2

    # create 1d arays of the possible x and y values (for later use)
    xs = np.arange(nx)
    ys = np.arange(ny)
    
    # pixel coordinates for each pixel (for later use)
    yvals, xvals = np.indices(image.shape)
    # compute the row averages and normalize so that the background is near 0 and the peaks are near 1
    rowaverage = image.mean(axis=1)
    rowaverage -= np.median(rowaverage)
    rowaverage /= rowaverage.max()
    # plt.plot(ys, rowaverage)
    # plt.xlabel('Row Number (y-coordinate)'), plt.ylabel('Normalized Row Average')
    # plt.grid()
    # plt.show()
    
    # #### Notes
    # * the plot above should show peaks in the rows containing object light
    # * you may also notice a few bad rows (actually bad columns!) with counts significantly below the median.
    # Based on the plot above, record the row coordinates (`ys`) that are brighter that some threshold, say 20% of the profile peak. To be conservative, also include the 5 rows above and 5 rows below.
    # 
    # Then create a 2D boolean mask that is `True` for pixels that are dominated by sky light and `False` for other pixels.
    
    # find the rows with object light
    objrows = ys[rowaverage > 0.2]
    
    # add some margin to object rows
    ngrow = 1 # number of rows to include above and below object rows
    newobjrows = []
    for row in objrows:
        newobjrows.extend([row + i for i in np.arange(-ngrow, ngrow + 1)])
    objrows = np.unique(newobjrows)[:-1]
    
    # mask to mark sky rows
    skymask = np.ones(image.shape, dtype=bool)
    print(objrows.shape,skymask.shape)
    # show_image(skymask)
    # plt.show()
    # input()
    skymask[objrows, :] = False
    
    # also exclude bad rows
    badrows = ys[rowaverage < -0.05]
    skymask[badrows, :] = False
    
    # rows with mostly sky background light
    skyrows = ys[skymask.mean(axis=1) == 1]
    
    # With the object traces and bad rows masked, we can also check for cosmic rays and mask these as well. To do this we calculate the median value and standard deviation along each column and then reject pixels that are more than a certain number (typically 3-5) of standard deviations away from the median. These deviant pixels are then noted on the `skymask`.
    # median (unmasked) sky spectrum and standard deviation
    medspec = np.median(image[skyrows, :], axis=0)
    stdspec = np.std(image[skyrows, :], axis=0, ddof=1)
    
    # exclude deviant pixels from the skymask
    pull = (image - medspec) / stdspec
    w = pull > 5
    skymask[w] = False
    
    # show the mask
    # plt.figure(figsize=(15, 3))
    # plt.imshow(skymask, origin='lower', aspect='auto')
    # plt.show()
    
    # ## Look at a small section of image up close

    # cut out a small image "stamp" near the center of the frame
    row = cy -400
    col = cx -400
    hwidth = 20
    
    image_raw = image
    image = ma.masked_array(image,np.logical_not(skymask))
    
    stamp = image[:, col - hwidth : col + hwidth]
    ys_stamp = yvals[:, col - hwidth : col + hwidth]
    xs_stamp = xvals[:, col - hwidth : col + hwidth]
    
    # show the image stamp
    extent = (xs_stamp.min(), xs_stamp.max(), ys_stamp.min(), ys_stamp.max())
    # show_image(stamp, extent=extent)
    # plt.show()
    
    # Recall the vertical bands are sky lines that mark lines of constant wavelength. Notice that these do not run perfectly parallel to the columns. Rather, the `x` coordinate for a given wavelength will be a function of the row number. 
    # 
    # As the plot should demonstrate, the lines of constant wavelength are slightly tilted with respect to the columns and there is also slight curvature. Thus we can approximate that the `x` coordinate for a given wavelength will be a quadratic function of the row number.
    # #### Note
    # Because wavelength varies along a given column, if we simply plot the counts in each pixel against each pixel's column number then we get a range of values in each column:

    # plot stamp values against column numbers
    # plt.plot(xs_stamp.ravel(), stamp.ravel(), 'r.');
    # plt.xlabel('Column Number'), plt.ylabel('Counts');
    # plt.show()
    
    # ## Map out lines of constant wavelength (determine the wavelength for each pixel)
    # We can model the change in wavelength coordinate in arbitrary units, `dl`, from the wavelength at some reference pixel in terms of the offsets, `dx` and `dy` from the reference pixel.
    # 
    # We can then write down a function that takes our model parameters (the slope and curvature of the lines of constant wavelength with respect to the columns), the offsets in the x and y coordinates, `dxs` and `dys`, respectively, and returns the wavelength offset from the reference pixel (i.e. `dx = dy = 0`).
    
    def get_dl_model(params, dxs, dys):
        return dxs + params[0] * dys + params[1] * dys ** 2
    
    # To see how this works, make a guess for the slope and curvature parameters and try plotting the wavelength offsets from our reference pixel (the one at `x=col` and `y=row`).
    # pixel offsets from the refernece pixel
    dxs = xs_stamp - col
    dys = ys_stamp - row
    
    # parameter guess
    guess = (-0.01, 0)
    
    # get the wavelength offsets and plot vs. counts
    dls = get_dl_model(guess, dxs, dys)
    # plt.plot(dls.ravel(), stamp.ravel(), 'r.')
    # plt.xlabel('Wavelength Offset')
    # plt.ylabel('Counts')
    # plt.show()
    
    
    # You should notice that for the right choice of model parameters, the vertical scatter is significantly reduced. This demonstrates one way to decide on the best model parameters: find the model parameters that minimize the vertical scatter. But how do we measure this scatter?
    # ## Fit a spline to the counts spectrum above
    # Given the data above (the wavelength offsets, `dls`, and the stamp count values) we can fit a spline to the resulting curve above and then subtract this off from the curve and measure the scatter. Notice that for every column we get a lot of samples (`ny` to be precise). Thus we can set the spline knots every pixel (or even closer). 
    # We can define a function that will return the best fit spline object (in the least-squares sense) for our data.

    def get_profile_spl(dls, stamp):
        # need to sort the data (and weights) so that the x values are increasing
        x, y = dls.ravel(), stamp.ravel()
        weights = np.sqrt(np.abs(y)) # not technically optimal for coadded data, but ok
        wsort = x.argsort()
        x, y, weights = x[wsort], y[wsort], weights[wsort]
    
        # set locations for spline knots
        t = np.linspace(x.min() + 1, x.max() - 1, np.int(x.max() - x.min()))
        spl = LSQUnivariateSpline(x, y, t, weights)
        return x, y, spl
    
    # fit a spline to the data and plot
    x, y, spl = get_profile_spl(dls, stamp)
    
    # fig, axarr = plt.subplots(2, sharex=True)
    # axarr[0].plot(x, y, 'r.')
    # axarr[0].plot(x, spl(x))
    # axarr[1].plot(x, y - spl(x), 'r.')
    # plt.ylim(-200, 200)
    # plt.xlabel('Wavelength Offset')
    # plt.show()
    
    # #### Notice
    #  * the counts may vary along the slit (both because of real spatial variations and, more importantly, we have not accounted for efficiency variations along the slit), so there may be some biased scatter
    #  * the parameter guess above may be imperfect resulting in systematic residuals (we can minimize these below)
    # ## Minimize residuals from the spline model to determine the best tilt/curvature model
    # We need a metric to determine how well the spline fits the data. We could take a simple sum of the squares of the residuals or weight these by the expected poison noise to determine a $\chi^2$ value.
    
    def check_dl_model(params, dxs, dys, stamp):
        dls = get_dl_model(params, dxs, dys)
        x, y, spl = get_profile_spl(dls, stamp)
        chisq = np.sum((stamp - spl(dls)) ** 2 / np.sqrt(np.abs(stamp)))
        return chisq / (stamp.size - len(params))
    
    
    # In[841]:
    
    
    # see how good our guess is
    check_dl_model(guess, dxs, dys, stamp)
    
    
    # Now we just need to change the model parameters to minimize the residuals. We can use `fmin` for this.
    
    # In[842]:
    
    
    # get the best model parameters for this stamp
    params = fmin(check_dl_model, guess, args=(dxs, dys, stamp))
    print("best model parameters are", params)
    
    
    # You can plug in these parameters for your `guess` above and re-run the cells to check the residuals if you want.
    
    # ## Find the lines of constant wavelength at other parts of the 2-D image
    
    # The tilt/curvature between lines of constant wavelength and the columns can vary across the spectrum. So far we have only determined it for a small portion of the image. Now we can apply the same techniques to map out the relative wavelengths shifts across the image.
    # 
    # To do this, we need to pick a width (number of columns) within which we can expect to get enough sky lines to do the mapping. Considering the left side (shorter wavelength or "blue" portion) of the spectrum, we can set the width at about 400 columns. (Note it is possible to vary the width across the image to take advantage of the extra wavelength information in the near-IR, but this is not required).
    
    # In[843]:
    
    
    # define the column centers for the stamps
    hwidth = 200
    cols = np.arange(hwidth, nx, 2 * hwidth)
    cols
    
    
    # As above, map out the wavelength shits along columns, one stamp (image section) at a time. It is best to use the sky mask we defined above and possibly some outlier rejection to fend off deviant pixels.
    # 
    # Notice that the code cell below will keep track of the wavelength offsets for all columns, but really we only need the central columns defined above.

    # reference wavelength offsets to the center row
    row = cy
    
    # define a 2D array to hold the wavelength offsets for each pixel
    lambdas = np.zeros(image.shape) 
    
    # loop over each central column
    for col in cols:
        print('col = ', col)
        
        # slice the data
        inds = np.s_[:, col - hwidth : col + hwidth]
        stamp = image[inds]
        mask = skymask[inds]
        dys = yvals[inds] - row
        dxs = xvals[inds] - col
    
        # initial fit
        params = fmin(check_dl_model, guess, args=(dxs[mask], dys[mask], stamp[mask]))
        
        # check for outliers
        dls = get_dl_model(guess, dxs, dys)
        x, y, spl = get_profile_spl(dls, stamp)
        model = spl(dls)
        pull = (stamp - model) / np.sqrt(np.abs(model))
        w = (pull < 5) & mask
        params2 = fmin(check_dl_model, params, args=(dxs[w], dys[w], stamp[w]))
    
        # record
        lambdas[inds] = get_dl_model(params2, dxs, dys) + col
    
    
    # Look at a few rows and see how the wavelength offsets vary with column number. We can fit a low-order polynomial to these.
    
    # In[845]:
    
    
    # just plot offsets for a few of the rows across the image
    order = 3
    # print (ny)
    for y in range(10, ny, 15):
        # p = plt.plot(cols, lambdas[y, cols] - xs[cols], 'o')
        c = np.polyfit(cols, lambdas[y, cols] - xs[cols], order)
    #     plt.plot(xs, np.polyval(c, xs), c=p[0].get_color(), label='row {}'.format(y))
    # plt.legend()
    # plt.xlabel('Column Number')
    # plt.ylabel('Wavelength Offset from Middle Row')
    # plt.show() 
    
    # You may notice that the tilt (the wavelength difference in a given column between the first row and the last row) increases as we move to larger column numbers. Make sure the order is large enough to follow the trends without over fitting the data (i.e. there should not be noticeable wiggles between data points).
    # 
    # Now we can fit for every row. We could do a 2D fit to all the data at once, but simply fitting row by row works here.
    
    # get the lambda values for the entire image (fit)
    lambdafit = np.zeros(image.shape)
    for y in range(ny):
        c = np.polyfit(cols, lambdas[y, cols] - xs[cols], order)
        lambdafit[y, :] = np.polyval(c, xs) + xs
        
    print (lambdafit.shape)


    # function to fit a 2D spline
    def fit_sky(xvals, yvals, image, ky=1, dx=0.5):
        # select knot points in the x (wavelength) direction
        tx = np.arange(xvals.min() + 2, xvals.max() - 2, dx)
        
        # select knot points in the y (spatial) direction
        ty = [] # if there are no knots, the fit will be a poly nomial of degree ky
        
        # fit the 2D spline
        return LSQBivariateSpline(xvals.ravel(), yvals.ravel(), image.ravel(), tx, ty, ky=ky)

     # use the (unmasked) sky background pixels and fit the 2D spline
    skyfit = fit_sky(lambdafit[skymask], yvals[skymask], image[skymask])
    
    # evaluate the 2D sky at every pixel
    sky = skyfit.ev(lambdafit, yvals)
    
    # print (sky.shape)

    # show_image(image)
    # plt.show()
    show_image(sky)
    plt.savefig(target_name+'_'+image_name.split('.')[0]+'_fitted_sky.png')
    # plt.show()
    # # plt.show()
    # show_image(image-sky)
    # plt.show()
    # input()

    sky_subtracted = np.array(image)-np.array(sky)
    # show_image(sky_subtracted)
    # plt.show()

    plt.close()
    fits.writeto(filename='s'+image_name,data=np.array(sky_subtracted),header=header,overwrite=True)

    return image_raw, sky_subtracted, sky, xs, ys, nx, ny, yvals


def optimal(image_raw, sky, xs, ys, nx, ny, yvals,trace_c,display):    

    trace_c = np.poly1d(trace_c[0]) 

    # print (trace_c)
    # input()
    def get_profile_model(params, ys):
        a1, cy1, sigma1 = params
        
    #     p1 = np.exp(-(ys - cy1)**2 / 2 / sigma1**2) 
    #     p1 /= p1.max()
        
    #     # p2 = np.exp(-(ys - cy2)**2 / 2 / sigma2**2) 
    #     # p2 /= p2.max()

    #     return a1 * p1 #+ a2 * p2

    # # get the median for each row
    # profile = np.median(image_raw - sky, axis=1)
    
    # # starting guess for the profile model
    # guess = (80, 31, 2)
    # model = get_profile_model(guess, ys)
    
    # plt.plot(ys, profile,label='Data')
    # plt.plot(ys, model,label='Model')
    # plt.xlabel('Row Number')
    # plt.ylabel('Median Row Counts')
    # plt.legend(loc='best')
    # plt.show()

    def get_profile_chisq(params, ys, profile):
        model = get_profile_model(params, ys)
        return np.sum( (profile - model)**2 / np.sqrt(np.abs(profile)) ) / (profile.size - len(params))
    # fit for the best model
    # params = fmin(get_profile_chisq, guess, args=(ys, profile))
    # print("best fit parameters are", params)
    
    # model = get_profile_model(params, ys)
    # plt.plot(ys, profile,label='Data')
    # plt.plot(ys, model,label='Model')
    # plt.xlabel('Row Number')
    # plt.ylabel('Median Row Counts')
    # plt.legend(loc='best')
    # plt.show()
    # input()

    # fit the profile centered at these columns
    hwidth = 50
    cols = np.arange(hwidth, nx + 1, 2 * hwidth)
    
    ycenter = np.zeros( (len(cols), 2) )
    # for icol, col in enumerate(cols):
    #     stamp = (image_raw - sky)[:, col - hwidth : col + hwidth]
    #     profile = np.mean(stamp, axis=1)
    #     params = fmin(get_profile_chisq, guess, args=(ys, profile))
    #     ycenter[icol, :] = params[[1, 2]]

    # fit the relation with a polynomial
    # ind = 0 # which trace 0 or 1?
    # t_order = 3
    # trace_c = np.polyfit(cols, ycenter[:, ind], t_order)
    # input()
    # fig, axarr = plt.subplots(2, sharex=True)
    # # axarr[0].plot(cols, ycenter[:, ind], 'ro')
    # axarr[0].plot(xs, np.polyval(trace_c, xs), 'r')
    # axarr[0].axes.set_ylabel('y-coordinate'); axarr[0].grid();
    # axarr[1].plot(cols, ycenter[:, ind] - np.polyval(trace_c, cols), 'ro')
    # axarr[1].axes.set_ylim(-0.5, 0.5)
    # axarr[1].axes.set_ylabel('Fit Residual (pixels)')
    # plt.xlabel('Column Number'); axarr[1].grid()
    # plt.show()

    # position offsets from the object trace (defined to be at slitpos = 0)
    slitpos = yvals - np.polyval(trace_c, yvals)

    # subtract the sky
    nosky = image_raw - sky
    
    # show_image(nosky)
    # plt.plot(cols, ycenter[:, ind], 'ro')
    # plt.plot(xs, np.polyval(trace_c, xs), 'r')
    # plt.show()
    
    
    # normalize to the pixel brightness at the trace center
    yinds = (np.round(np.polyval(trace_c, xs))).astype(int)
    normed = nosky / nosky[yinds, xs]

    # get 1D arrays with the positions along the slit and the normalized counts
    pos = slitpos.flatten()
    counts = normed.flatten()
    
    # sort by slit position
    sort_inds = pos.argsort()
    pos, counts = pos[sort_inds], counts[sort_inds]
    
    # fit a spline to model the spatial profile
    t = np.linspace(pos.min() + 2, pos.max() - 2, ny // 2) # spline knot points
    profile_spl = LSQUnivariateSpline(pos, counts, t)
    
    # remove outliers and re-fit
    diff = counts - profile_spl(pos)
    sample = sigma_clip(diff)
    w = ((np.abs(diff) / sample.std()) < 5) & np.isfinite(diff)
    profile_spl = LSQUnivariateSpline(pos[w], counts[w], t)

    # plot the target profile
    # plt.plot(pos, profile_spl(pos) )
    # plt.xlim(-40, 50)
    # plt.grid()
    # plt.show()

    # create the profile image
    profile_image = profile_spl(slitpos)
    
    # de-weight negative values in provile_image
    profile_image[profile_image < 0] = 0

    
    if display is True:
        show_image(profile_image, upper=50)
        plt.show()

    # select which rows to sum
    w = (slitpos > -10) & (slitpos < 10)
    ymin, ymax = yvals[w].min(), yvals[w].max()
    
    print (slitpos.shape)
    
    # calculate the sum
    spec_basic = nosky[ymin:ymax, :].sum(axis=0)
    
    # sky background
    skybg_basic = sky[ymin:ymax, :].sum(axis=0)
    
    # plot the extracted spectrum
    # plt.plot(xs, spec_basic)
    # plt.xlabel('Column Number')
    # plt.ylabel('Counts')
    # plt.show()

    # calculate the weighted average (for each column)
    spec_opt = (nosky * profile_image)[ymin:ymax, :].sum(axis=0) / profile_image.sum(axis=0)
    
    # calculate the bias factor needed to scale the average to a sum
    bias_factor = np.median(spec_basic / spec_opt)
    spec_opt *= bias_factor
    
    # same for the sky background
    skybg_opt = (sky * profile_image)[ymin:ymax, :].sum(axis=0) / profile_image.sum(axis=0)
    bias_factor_sky = np.median(skybg_basic / skybg_opt)
    skybg_opt *= bias_factor_sky
    
    # plot the extracted spectrum
    # plt.plot(xs, spec_basic, label='basic extraction')
    # plt.plot(xs, spec_opt, label='optimal extraction')
    # plt.xlabel('Column Number')
    # plt.ylabel('Counts')
    # plt.show()

    return spec_opt



def ap_trace(image, object_keyword,
    trace_prominence, y0_trace, yf_trace, 
    manual_trace, manual_poly_order, manual_x, manual_y, 
    tolerance, trace_width, fmask=(1,), nsteps=20,
    recenter=False, prevtrace=(0,), bigbox=15, Saxis=1, poly_order=2, 
    display=False):
    """
    Trace the spectrum aperture in an image
    Assumes wavelength axis is along the X, spatial axis along the Y.
    Chops image up in bins along the wavelength direction, fits a Gaussian
    within each bin to determine the spatial center of the trace. Finally,
    draws a cubic spline through the bins to up-sample the trace.
    Parameters
    ----------
    image : 2d numpy array
        This is the image, stored as a normal numpy array. Can be read in
        using astropy.io.fits like so:
        >>> hdu = fits.open('file.fits')  # doctest: +SKIP
        >>> img = hdu[0].data  # doctest: +SKIP
    object_keyword: Keyword defined in header that indicates target's name
    trace_prominence: float 
        Sets level to search for trace peak.  
    y0_trace : pixel number (int),
        The centre of the trace (y-value) at x = 0.  Sets start of trace position.
        Especially important when there is more than one trace on the ccd.
    yf_trace : pixel number (int),
        The centre of the trace (y-value) at the end of the trace.  Sets end of trace position.
        Especially important when there is more than one trace on the ccd.
    tolerance : int, default = 3
        Sets the dy limit within which to search for trace.
    trace_width :  float (pixel space)
        If trace width is None, then default is +- fwhm of trace fit. 
    Poly_order : Integer
        Order of polynomial that will be fitted to trace. Default is 2
    nsteps : int, optional
        Keyword, number of bins in X direction to chop image into. Use
        fewer bins if ap_trace is having difficulty, such as with faint
        targets (default is 50, minimum is 4)
    fmask : array-like, optional
        A list of illuminated rows in the spatial direction (Y), as
        returned by flatcombine.
    interac : bool, optional
        Set to True to have user click on the y-coord peak. (Default is
        False)
    recenter : bool, optional
        Set to True to use previous trace, but allow small shift in
        position. Currently only allows linear shift (Default is False)
    bigbox : float, optional
        The number of sigma away from the main aperture to allow to trace
    display : bool, optional
        If set to true display the trace over-plotted on the image
    Saxis : int, optional
        Set which axis the spatial dimension is along. 1 = Y axis, 0 = X.
        (Default is 1)
    Returns
    -------
    my : array
        The spatial (Y) positions of the trace, interpolated over the
        entire wavelength (X) axis
    """
    cwd = os.getcwd()
    data_dir = cwd

    # if subtract_sky == False:
    header = fits.getheader(image)
    img = fits.getdata(image)
    # img, header = fits.getdata(os.path.join(data_dir,image), header=True)   

    target_name = header[object_keyword]
    file_name = image.split('.')[0]

    if display is True:
        show_image(img)
        plt.title('Trace region of ' + target_name)
        # plt.axhline(trace_ymin)
        # plt.axhline(trace_ymax)
        plt.savefig(target_name+'_'+ file_name +'_trace_region.png')
        plt.show()

    if display is False:
        plt.title('Trace region of ' + target_name)
        show_image(img)
        # plt.axhline(trace_ymin)
        # plt.axhline(trace_ymax)
        plt.savefig(target_name+'_'+ file_name +'_trace_region.png')
        plt.close()


    img_window = img[:,:]
    # img_window = img[trace_ymin:trace_ymax,:]

    # show_image(img_window)
    # plt.show()

    # define the wavelength axis
    Waxis = 0
    # add a switch in case the spatial/wavelength axis is swapped
    if Saxis == 0:
        Waxis = 1

    print('Tracing Aperture using nsteps='+str(nsteps))
    # the valid y-range of the chip
    if (len(fmask)>1):
        ydata = np.arange(img_window.shape[Waxis])[fmask]
    else:
        ydata = np.arange(img_window.shape[Waxis])

    # need at least 4 samples along the trace. sometimes can get away with very few
    if (nsteps<4):
        nsteps = 4

    # median smooth to crudely remove cosmic rays
    img_sm = scipy.signal.medfilt2d(img_window, kernel_size=(5,5))
    # img_sm = img

    #--- Pick the strongest source, good if only 1 obj on slit
    ztot = img_window.sum(axis=Saxis)[ydata]
    yi = np.arange(img_window.shape[Waxis])[ydata]
    peak_y = yi[np.nanargmax(ztot)]
    # peak_y = 7.5
    peak_guess = [np.nanmax(ztot), np.nanmedian(ztot), peak_y, 2.]

    #-- use middle of previous trace as starting guess
    if (recenter is True) and (len(prevtrace)>10):
        peak_guess[2] = np.nanmedian(prevtrace)

    peaks, _ = find_peaks(ztot[np.isfinite(ztot)], prominence=trace_prominence)
    # plt.plot(yi[np.isfinite(ztot)], ztot[np.isfinite(ztot)])
    # plt.plot(peaks, ztot[np.isfinite(ztot)][peaks], "x")
    # plt.show()
    # input()

    popt_tot, pcov = curve_fit(_gaus, yi[np.isfinite(ztot)], ztot[np.isfinite(ztot)], p0=peak_guess)
    
    #-- only allow data within a box around this peak
    ydata2 = ydata[np.where((ydata>=popt_tot[2] - popt_tot[3]*bigbox) &
                            (ydata<=popt_tot[2] + popt_tot[3]*bigbox))]

    ztot = img_window.sum(axis=Saxis)[ydata2]
    yi = np.arange(img_window.shape[Waxis])[ydata2]
    # define the X-bin edges
    xbins = np.linspace(0, img_window.shape[Saxis], nsteps, dtype='int')
    ybins = np.zeros_like(xbins, dtype='float')
    fwhm = np.zeros_like(xbins, dtype='float')

    # Manual trace

    if manual_trace is True:
        def man_fit(manual_x,manual_y,manual_poly_order):
            ap_spl_man = np.polyfit(manual_x, manual_y, deg=manual_poly_order,full=True)
            p_man = np.poly1d(ap_spl_man[0])
            my_man = p_man(xbins) #+ trace_ymin
            show_image(img)
            plt.plot(xbins,my_man,'--r')
            plt.plot(manual_x,manual_y,'or')
            plt.show()
            return ap_spl_man, p_man, my_man
        ap_spl_man, p_man, my_man = man_fit(manual_x,manual_y,manual_poly_order)
        adjust = input('Would you like to adjust trace positions? (y/n) ')
        while adjust == 'y':
            manual_x = [float(kk) for kk in input('Enter x values: (eg. 0 1000 ... 1700) ').split()]
            manual_y = [float(kk) for kk in input('Enter y values: (eg. 25 27 ... 30) ').split()]
            ap_spl_man, p_man, my_man = man_fit(manual_x,manual_y,manual_poly_order)
            adjust = input('Would you like to adjust trace positions? (y/n) ')
        

    def straight(x,m,c):
        y = m*x + c
        return y

    slope, intercept, r_value, p_value, std_err = linregress([xbins[0],xbins[-1]],[y0_trace,yf_trace])
    yvals_line = straight(xbins,slope,intercept)
    
    plt.figure(figsize=(10, 7))
    for i in range(0,len(xbins)-1):
        #-- fit gaussian w/i each window
        if Saxis == 1:
            zi = np.sum(img_window[ydata2, xbins[i]:xbins[i+1]], axis=Saxis)
        else:
            zi = img_window[xbins[i]:xbins[i+1], ydata2].sum(axis=Saxis)

        # plt.clf()
        peaks, _ = find_peaks(zi[np.isfinite(ztot)], prominence=trace_prominence,threshold=tolerance)
        # plt.plot(yi[np.isfinite(ztot)],zi[np.isfinite(ztot)])
        # plt.plot(yi[np.isfinite(ztot)][peaks], zi[np.isfinite(ztot)][peaks], "x")
        # plt.plot(yi[np.isfinite(ztot)],zi[np.isfinite(ztot)],'r')
        # plt.title(xbins[i])
        # plt.show()

        # input('Enter')

        trace_start = min(yi[np.isfinite(ztot)][peaks], key=lambda x:abs(x-yvals_line[i])) 
        # tolerance = 2
    
        if i == 0:
            # inten = np.nansum(img_window[xbins[i],trace_start-tolerance:trace_start+tolerance])
            pguess = [np.abs(np.nanmax(zi)),np.abs( np.nanmedian(zi)), trace_start, 2.]
            popt,pcov = curve_fit(_gaus, yi[np.isfinite(ztot)], zi[np.isfinite(ztot)], p0=pguess, 
                bounds=(0,[abs(np.nanmax(zi)*10), abs(np.nanmedian(zi)*10), trace_start+tolerance, 4.]))
            if abs(yvals_line[i]-popt[2]) > tolerance:
                y1_trace = yvals_line[i]
            else:
                y1_trace = popt[2]
            # print (trace_start,y1_trace,popt[2])
            # input()
        if i >= 1:
            trace_next = min(yi[np.isfinite(ztot)][peaks], key=lambda x:abs(x-yvals_line[i])) 
            pguess = [np.abs(np.nanmax(zi)),np.abs( np.nanmedian(zi)), trace_next, 2.]
            popt,pcov = curve_fit(_gaus, yi[np.isfinite(ztot)], zi[np.isfinite(ztot)], p0=pguess,
                bounds=(0,[abs(np.nanmax(zi)*10), abs(np.nanmedian(zi)*10), trace_start+tolerance, 4.]))
                #  bounds=(([abs(np.nanmax(zi)*.1), abs(np.nanmedian(zi)*.1), 0, 0]),
                # ([abs(np.nanmax(zi)*10), abs(np.nanmedian(zi)*10), trace_next+tolerance, 4.])))
            # print (trace_next,popt[2])
            if abs(yvals_line[i]-popt[2]) > tolerance:
                y1_trace = yvals_line[i]
                # y0_trace = yvals_line[i]
            else:
                y1_trace = popt[2]
                # y0_trace = popt[2]
        
        plt.subplot(1, 2, 1)
        plt.plot(yi[np.isfinite(ztot)],zi[np.isfinite(ztot)]) # ,label=str(xbins[i])+' - '+str(np.round(popt[2],2))
        if i == int(nsteps*.5):
            plt.axvline(popt[2]-tolerance,color='k',ls='--')
            plt.axvline(popt[2]+tolerance,color='k',ls='--')
        plt.xlabel('Row Number')
        plt.ylabel('Intensity')
        
        # plt.plot(yi[np.isfinite(ztot)],_gaus(yi[np.isfinite(ztot)],*popt),'ro:',label='fit')
        
        ybins[i] = y1_trace #popt[2]
        fwhm[i] = popt[3]*gaussian_sigma_to_fwhm 

    # recenter the bin positions, trim the unused bin off in Y

    # input()

    Mxbins = (xbins[:-1]+xbins[1:]) / 2.
    Mybins = ybins[:-1]

    if trace_width == None:
        myfwhm = max(fwhm)*1.2
        # print (fwhm, myfwhm)
        # input()
    else:
        myfwhm = trace_width

    # mxbins[0] = 0
    # mybins[0] =  5.29467085
    # print (mxbins,mybins)

    
    def distance(x,y,tolerance):
        mxbins = []
        mybins = []
        for i in range(len(x)-1):
            d = y[i+1] - y[i] 
            # print (d)
            if d < .5*tolerance:
                mxbins.append(x[i])
                mybins.append(y[i])
        return mxbins, mybins

    mxbins, mybins = distance(Mxbins,Mybins,tolerance)


    ap_spl = np.polyfit(mxbins, mybins, deg=poly_order,full=True)
    print ('residuals = ',ap_spl[4])

    p = np.poly1d(ap_spl[0])

    # interpolate the spline to 1 position per column
    mx = np.arange(0, img_window.shape[Saxis])
    my = p(mx) #+ trace_ymin




    if display is True:
        plt.title('Gaussian fitted to trace of ' + target_name)
        # plt.legend(loc='best')
        plt.subplot(1, 2, 2)
        plt.title('Fitted centers of trace')
        plt.plot(mxbins,mybins,'k.')
        plt.plot(mxbins,p(mxbins),'r:',label='fit')
        res = '%.1e' % ap_spl[4]
        plt.plot(mxbins,p(mxbins),'r:',label='residuals = '+ str(res))
        plt.xlabel('Column Number')
        plt.ylabel('Row Number')
        plt.legend(loc='best')
        plt.savefig(target_name+'_'+ file_name +'_trace_fit.png')
        plt.show()

        show_image(img)
        plt.autoscale(False)
        plt.plot(mx,my,'b',lw=1,label='Fitted trace')
        plt.plot(mx,my - myfwhm,'k',label='Trace width')
        plt.plot(mx,my + myfwhm,'k')
        if manual_trace is True:
            plt.plot(xbins,my_man,'r',label='Manual trace')
            plt.plot(manual_x,manual_y,'or')
        plt.legend(loc='best')
        plt.title('Trace region of ' + target_name)
        plt.savefig(target_name+'_'+ file_name +'_trace_fit2D.png')
        plt.show()
    
    if display is False:
        plt.title('Gaussian fitted to trace of ' + target_name)
        # plt.legend(loc='best')
        plt.subplot(1, 2, 2)
        plt.title('Fitted centers of trace')
        plt.plot(mxbins,mybins,'k.')
        plt.plot(mxbins,p(mxbins),'r:',label='fit')
        res = '%.1e' % ap_spl[4]
        plt.plot(mxbins,p(mxbins),'r:',label='residuals = '+ str(res))
        plt.xlabel('Column Number')
        plt.ylabel('Row Number')
        plt.legend(loc='best')
        plt.savefig(target_name+'_'+ file_name +'_trace_fit.png')
        plt.close()

        show_image(img)
        plt.autoscale(False)
        plt.plot(mx,my,'b',lw=1,label='Fitted trace')
        plt.plot(mx,my - myfwhm,'k',label='Trace width')
        plt.plot(mx,my + myfwhm,'k')
        if manual_trace is True:
            plt.plot(xbins,my_man,'r',label='Manual trace')
            plt.plot(manual_x,manual_y,'or')
        plt.legend(loc='best')
        plt.title('Trace region of ' + target_name)
        plt.savefig(target_name+'_'+ file_name +'_trace_fit2D.png')
        plt.close()

    plt.close()
    # print("> Trace gaussian width = "+str(popt_tot[3])+' pixels')
    # print (my, myfwhm)
    if manual_trace is True:
        return my, myfwhm, ap_spl, my_man, ap_spl_man
    else:    
        return my, myfwhm, ap_spl

# my, myfwhm, trace_c = ap_trace(IMAGE, fmask=(1,), nsteps=5, interac=False,
#              recenter=True, prevtrace=(0,), bigbox=15,
#              Saxis=1, display=True,object_keyword='OBJECT')


def ap_extract(image, trace, object_keyword, gain_keyword, readnoise_keyword, 
    apwidth=8, skysep=3, skywidth=7, skydeg=0,coaddN=1,display=False):
    """
    1. Extract the spectrum using the trace. Simply add up all the flux
    around the aperture within a specified +/- width.
    Note: implicitly assumes wavelength axis is perfectly vertical within
    the trace. An major simplification at present. To be changed!
    2. Fits a polynomial to the sky at each column
    Note: implicitly assumes wavelength axis is perfectly vertical within
    the trace. An important simplification.
    3. Computes the uncertainty in each pixel
    Parameters
    ----------
    img : 2d numpy array
        This is the image, stored as a normal numpy array. Can be read in
        using astropy.io.fits like so:
        >>> hdu = fits.open('file.fits') # doctest: +SKIP
        >>> img = hdu[0].data # doctest: +SKIP
    trace : 1-d array
        The spatial positions (Y axis) corresponding to the center of the
        trace for every wavelength (X axis), as returned from ap_trace
    apwidth : int, optional
        The width along the Y axis on either side of the trace to extract.
        Note: a fixed width is used along the whole trace.
        (default is 8 pixels)
    skysep : int, optional
        The separation in pixels from the aperture to the sky window.
        (Default is 3)
    skywidth : int, optional
        The width in pixels of the sky windows on either side of the
        aperture. (Default is 7)
    skydeg : int, optional
        The polynomial order to fit between the sky windows.
        (Default is 0)
    Returns
    -------
    onedspec : 1-d array
        The summed flux at each column about the trace. Note: is not
        sky subtracted!
    skysubflux : 1-d array
        The integrated sky values along each column, suitable for
        subtracting from the output of ap_extract
    fluxerr : 1-d array
        the uncertainties of the flux values
    """

    cwd = os.getcwd()
    data_dir = cwd

    if image[0] == 's':
        sky = True
    else:
        sky = False

    file_name = image.split('.fits')[0]
    image, header = fits.getdata(os.path.join(data_dir,image), header=True) 
    target_name = header[object_keyword]

    onedspec = np.zeros_like(trace)
    variancespec = np.zeros_like(trace)
    skysubflux = np.zeros_like(trace)
    fluxerr = np.zeros_like(trace)

    for i in range(0,len(trace)):
        #-- first do the aperture flux
        # juuuust in case the trace gets too close to the edge
        widthup = apwidth
        widthdn = apwidth
        if (trace[i]+widthup > image.shape[0]):
            widthup = image.shape[0]-trace[i] - 1
        if (trace[i]-widthdn < 0):
            widthdn = trace[i] - 1

        # simply add up the total flux around the trace +/- width
        onedspec[i] = image[int(trace[i]-widthdn):int(trace[i]+widthup+1), i].sum()

        #-- now do the sky fit
        itrace = int(trace[i])
        y = np.append(np.arange(itrace-apwidth-skysep-skywidth, itrace-apwidth-skysep),
                      np.arange(itrace+apwidth+skysep+1, itrace+apwidth+skysep+skywidth+1))

        z = image[y,i]

        if sky is not True:

            if (skydeg>0):
                # fit a polynomial to the sky in this column
                pfit = np.polyfit(y,z,skydeg)
                # define the aperture in this column
                ap = np.arange(trace[i]-apwidth, trace[i]+apwidth+1)
                # evaluate the polynomial across the aperture, and sum
                skysubflux[i] = np.sum(np.polyval(pfit, ap))
            elif (skydeg==0):
                skysubflux[i] = np.nanmean(z)*(apwidth*2.0 + 1)


        #-- finally, compute the error in this pixel
        sigB = np.std(z) # stddev in the background data
        N_B = len(y) # number of bkgd pixels
        N_A = apwidth*2. + 1 # number of aperture pixels

        # based on aperture phot err description by F. Masci, Caltech:
        # http://wise2.ipac.caltech.edu/staff/fmasci/ApPhotUncert.pdf


        if sky is True:
            fluxerr[i] = np.sqrt(np.sum((onedspec[i])/coaddN) +
                                 (N_A + N_A**2. / N_B) * (sigB**2.))
        else:
            fluxerr[i] = np.sqrt(np.sum((onedspec[i]-skysubflux[i])/coaddN) +
                                 (N_A + N_A**2. / N_B) * (sigB**2.))

        #-- calculate variance spectrum
        if type(gain_keyword) == str:
            gain = float(header[gain_keyword])  #1.145
        elif type(gain_keyword) == float:
            gain = gain_keyword
        if type(readnoise_keyword) == str:
            readnoise = float(header[readnoise_keyword]) 
        elif type(readnoise_keyword) == float:
            readnoise = readnoise_keyword
        # gain = float(header[gain_keyword])  #1.145
        # readnoise = float(header[readnoise_keyword])  # 2.245
        variancespec[i] = readnoise + image[int(trace[i]-widthdn):int(trace[i]+widthup+1), i].sum()/gain 


    # print (test,test.shape)
        
    # plt.plot(onedspec,label='spec')
    # plt.plot(skysubflux,label='sky')

    if sky is True:
        smooth_spec = smooth(onedspec,3)
    else:
        smooth_spec = smooth(onedspec-skysubflux,3)

    if display is True:
        plt.figure(figsize=(10, 5))
        plt.subplot(1,2,1)
        if sky is False:
            plt.plot(skysubflux,label='sky')
        plt.plot(smooth_spec,label='spec')
        plt.legend(loc='best')
        plt.title('Extracted spectrum: '+ target_name)
        plt.xlabel('Columns')
        plt.ylabel('Intesity')
        plt.subplot(1,2,2)
        plt.title('S/N spectrum: '+ target_name)
        if sky is False:
            plt.plot((onedspec-skysubflux)/np.sqrt(variancespec))
        else:
            plt.plot((onedspec)/np.sqrt(variancespec))
        plt.ylabel('S/N')
        plt.xlabel('Columns')
        plt.savefig(target_name+'_'+file_name+'_Initial_1D_spectrum.png')
        plt.show()

    if display is False:
        plt.figure(figsize=(10, 5))
        plt.subplot(1,2,1)
        if sky is False:
            plt.plot(skysubflux,label='sky')
        plt.plot(smooth_spec,label='spec')
        plt.legend(loc='best')
        plt.title('Extracted spectrum:'+ target_name)
        plt.xlabel('Columns')
        plt.ylabel('Intesity')
        plt.subplot(1,2,2)
        plt.title('S/N spectrum: '+ target_name)
        if sky is False:
            plt.plot((onedspec-skysubflux)/np.sqrt(variancespec))
        else:
            plt.plot((onedspec)/np.sqrt(variancespec))
        plt.ylabel('S/N')
        plt.xlabel('Columns')
        plt.savefig(target_name+'_'+file_name+'_Initial_1D_spectrum.png')
        plt.close()

    plt.close()
    if sky is True:
        return onedspec, fluxerr, variancespec
    if sky is False:
        return onedspec-skysubflux, fluxerr, variancespec


def redshift(w,z= 0.00093):    
    w_cor = w*(1+z)
    return w_cor


def wavelength(onedspec_sum,onedspec_optimal,spec_file_name,wave_min,wave_max,arc_file,reference_spec,line_list,trace,trace_fwhm,
    prominence,parameter_file,order,object_keyword,flip_wave_axis=False,view_arc=False,display=False):
    """
    1. Apply wavelength solution to onedspec from apextract.
    2. Fits a polynomial to interpolate between supplied column numbers and wavelength space
 
    Parameters
    ----------
    onedspec_sum : 1D numpy array
        This is the output from apextract for sumation across trace
    onedspec_optimal : 1D numpy array
        This is the output from optimal for optimal extraction
    wave_min : float
        Lower limit of wavelength range in Angstroms
    wave_max : float
        Upper limit of wavelength range in Angstroms
    arc_file : fits arc spectrum
    reference_spec : reference spectrum
        Data file with wavelength column in Angstroms (first column) and counts column (second column). 
        The file should not have any column headers.
    line_list : laboratory line list
        File with laboratory line list in Angstroms.  
        The first column contain wavelengts and second column contains line ID's (e.g CuI). 
        The file should not have any column headers.
    trace : 1D numpy array
        Trace returned by ap_trace.  
        Used to extract arc spectrum of same region as trace.
    trace_fwhm : float/int
        Average FWHM of fitted trace from ap_trace.  
        Used to extract arc spectrum of same region as trace.
    prominence : int
        Used by find_peaks (python function) to find spectral peaks
    parameter_file:  File 
        File with key parameters e.g. identified columns and wavlengths
        for wavelength solution.  See example file.  Best to keep format
        exactly as example.
    order :  int
        Polynomial fitted to interpolate between identified columns 
        and wavelenths.
    object_keyword :  Header keyword identifier
        Telescope specific, the header keyword that lists the object name.
    view_arc : True / False
        True: When you would like to look at the arc spectrum to idnetify column 
        numbers and wavelengths interactively. 
        False: Provide column numbers and corresponding wavelengths in parameter file.
    display : True / False
        If true plots will appear and be saved, if False plots will only be saved.

    Returns
    -------
    wavelength calibrated spectrum : two 1D arrays

    """

    global interact

    interact = view_arc

    spec_sum = np.flip(onedspec_sum)
    if onedspec_optimal is not None:
        spec_optimal = np.flip(onedspec_optimal)

    cwd = os.getcwd()
    data_dir = cwd
    row1 = int(np.median(trace)-trace_fwhm)
    row2 = int(np.median(trace)+trace_fwhm)

    file_name = arc_file.split('.fits')[0]
    header = fits.getheader(os.path.join(data_dir,arc_file))
    lamp_image = fits.getdata(os.path.join(data_dir,arc_file))
    # lamp_image, header = fits.getdata(os.path.join(data_dir,arc_file), header=True) 
    if flip_wave_axis is True:
        lamp_spec = np.flip(np.median(lamp_image[row1:row2, :], axis=0))  # 6965 A - closest to Halpha
    else:
        lamp_spec = np.median(lamp_image[row1:row2, :], axis=0)

    # plt.plot(np.arange(len(lamp_spec)),lamp_spec)
    # plt.show()

    peaks, _ = find_peaks(lamp_spec,prominence=prominence)

    target_name = image_header(spec_file_name+'.fits',object_keyword)

    # load the "linelist" with precise laboratory wavelengths for these lamp lines
    dtype = [('wav', float), ('id', 'U2')]
    linelist = np.genfromtxt(os.path.join(data_dir, line_list), dtype=dtype)
    linelist.sort(order='wav')

    # use a Gaussian function to model spectral lines
    def gaussian(x, *params):
        amp, x0, sigma = params
        return amp * np.exp(-(x - x0)**2 / 2 / sigma**2)


    # define a function to determine the precise column centers for each lamp line
    def get_lamp_lines(lamp_spec, prominence):
        peak_cols, properties = find_peaks(lamp_spec, prominence=prominence)
    
        # record in a structured array for convenience 
        dtype = []
        dtype.append( ('col', float) )
        dtype.append( ('counts', float) )
        dtype.append( ('x', float) )
        dtype.append( ('y', float) )
        dtype.append( ('sigma', float) )
        dtype.append( ('id', 'U2') )
        dtype.append( ('wav', float) )
        dtype.append( ('wavres', float) )
        dtype.append( ('used', bool) )
        lamp_lines = np.zeros(peak_cols.size, dtype=dtype)
        lamp_lines['col'] = peak_cols
        lamp_lines['counts'] = lamp_spec[peak_cols]
        lamp_lines['x'] = np.nan
        lamp_lines['y'] = np.nan
        lamp_lines['sigma'] = np.nan
        lamp_lines['wav'] = np.nan
        lamp_lines['wavres'] = np.nan
    
        # fit each peak to determine precise center
        cols = np.arange(lamp_spec.size)
        sigma_guess = 2.5
        for line in lamp_lines:
            if line['counts'] > 60000:
                # line is saturated...skip
                continue
    
            i0 = max([0, int(line['col'] - 5)])
            i1 = min([lamp_spec.size - 1, int(line['col'] + 5)])
            guess = (line['counts'], line['col'], sigma_guess)
            bounds = ((0, line['col'] - 3, 0), (np.inf, line['col'] + 3, np.inf))
            try:
                popt, pcov = curve_fit(gaussian, cols[i0:i1], lamp_spec[i0:i1], p0=guess, bounds=bounds)
            except RuntimeError:
                # curve_fit failed to converge...skip
                continue
    
            line['x'] = popt[1]
            line['y'] = gaussian(popt[1], *popt)
            line['sigma'] = popt[2]
    
        # filter lamp_lines to keep only lines that were fit
        wasfit = np.isfinite(lamp_lines['x'])
        lamp_lines = lamp_lines[wasfit]
        # print('found center pixel values for', lamp_lines.size, 'lines')
        return lamp_lines

    # helper function to mark lamp lines with and 'x'
    def mark_peaks(plt, lamp_lines, xtype='x', ytype='y', c='k'):
        w = np.isfinite(lamp_lines['wav'])
        if w.sum() > 0:
            plt.scatter(lamp_lines[xtype][w], lamp_lines[ytype][w], c=np.abs(lamp_lines['wavres'][w]), zorder=10)
        if (~w).sum() > 0:
            plt.scatter(lamp_lines[xtype][~w], lamp_lines[ytype][~w], c=c, marker='x')

    # find the precise column centers for all lamp lines
    lamp_lines = get_lamp_lines(lamp_spec, prominence=prominence)

    # print (lamp_lines)
    # input('Enter:')

    # load and plot the provided spectral atlas
    lamp_ref = np.genfromtxt(os.path.join(data_dir, reference_spec), names='wav, counts')

    # this function will be used below to match a given line to a given list
    def match_to_list(listing, values, plt=None, tol=None, revcoeff=None, c='k'):
        matched = []
        cols = []
        # print (listing)
        for value in values:
            # print (value)
            absdiff = np.abs(value - listing)
            ind = np.argmin(absdiff)
            if tol is None:
                bestmatch = listing[ind]
            elif absdiff[ind] < tol:
                bestmatch = listing[ind]
            else:
                bestmatch = np.nan
            matched.append(bestmatch)
    
            if plt is not None:
                plt.axvline(bestmatch, ls='dotted', c=c)
                
            if revcoeff is not None:
                col = np.polyval(revcoeff, bestmatch)
                cols.append(col)
                print(f"{bestmatch:.1f} is expected near column {col:.0f}")
    
        if revcoeff is not None:
            # print ('y',np.array(matched), cols)
            return np.array(matched), cols
        
        # print (matched)
        # input('Enter 0')
        return np.array(matched)

    wav1 = wave_min 
    wav2 = wave_max 

    # plot the lamp spectra
    if interact is True:
        rough_cols_old = []
        rough_waves_old = []

        def check_columns(rough_cols_old,rough_waves_old):
            plt.plot(lamp_spec, c='g')
            plt.xlabel('Column Number')
            plt.ylabel('Counts')
            plt.title('Note down columns of a few lines to be used for calibration:')
            plt.show()

            rough_waves = input('Enter rough wavelengths of spectral lines in Angstroms (e.g. 4198 4764 ... 6965): ').split(' ')
            rough_cols = input('Enter corresponding column numbers (e.g. 25 100 ... 2000): ').split(' ')

            rough_waves = [float(i) for i in rough_waves]
            rough_cols = [float(i) for i in rough_cols]

            plt.plot(lamp_spec, c='g')
            if len(rough_cols_old) != 0:
                for ii in range(len(rough_cols_old)):
                    plt.axvline(int(rough_cols_old[ii]),ls='--',color='r',lw=1.)
                    plt.text(int(rough_cols_old[ii])+10,max(lamp_spec),str(rough_waves_old[ii])+' $\AA$')
            plt.xlabel('Column Number')
            plt.ylabel('Counts')
            plt.title('Check identifications:')
            for ii in range(len(rough_cols)):
                plt.axvline(int(rough_cols[ii]),ls='--',color='k',lw=1.)
                plt.text(int(rough_cols[ii])+10,max(lamp_spec),str(rough_waves[ii])+' $\AA$')
            plt.show()
            recheck = input('Do you want to adjust columns and/or wavelengths? (y/n) ')
            return recheck, rough_waves, rough_cols

        recheck, rough_waves, rough_cols = check_columns(rough_cols_old,rough_waves_old)
        while recheck == 'y':
            rough_waves_old = rough_waves
            rough_cols_old = rough_cols
            for ii in range(len(rough_cols_old)):
                plt.axvline(int(rough_cols_old[ii]),ls='--',color='r',lw=1.)
                plt.text(int(rough_cols_old[ii])+10,max(lamp_spec),str(rough_waves_old[ii])+' $\AA$')
            recheck, rough_waves, rough_cols = check_columns(rough_cols_old,rough_waves_old)
    plt.close()
    if interact == True:
        inter = input('Do you want to continue interactively? (y/n) ')
        if inter == 'y':
            interact = True
        else:
            interact = False
            f = open(parameter_file,'w')
            f.write('rough_cols rough_waves\n')
            for ii in range(len(rough_cols)):
                f.write('%g %g\n' % (rough_cols[ii],rough_waves[ii]))
            f.close()



    if (interact == False) and parameter_file != None:
        f = open(parameter_file,'r')
        rough_cols, rough_waves = np.loadtxt(parameter_file,skiprows=1,usecols=(0,1),unpack=True)
        # for line in f:
        #     if line[0] != '#':
        #         data = list(line.split('=')[1].rstrip().split(' '))
        #         # print ('ok')
        #         if line.split(' = ')[0].rstrip() == 'rough_waves':
        #             print (line)
        #             rough_waves = [float(i) for i in data if len(i) != 0]

        #         if line.split('=')[0].rstrip() == 'rough_cols':
        #             rough_cols = [float(i) for i in data if len(i) != 0]
    
    refwavs = match_to_list(linelist['wav'], rough_waves, plt=plt)
    # print (linelist['wav'], refwavs)

    # input()

    # plt.show()

    # find a section of lamp_spec that looks similar
    # col1 = 10
    # col2 = 1900
    # plt.plot(lamp_spec, c='g')
    # plt.xlim(col1, col2)
    # plt.xlabel('Column Number'); plt.ylabel('Counts'); plt.grid()
    
    # # mark lines with Gaussian-fit centers
    # mark_peaks(plt, lamp_lines)
    # plt.show()

    # input('Enetr 2:')

    # record the rough column numbers of the same lines as above in the same order
    # rough_cols = [230,427,1194,1382,1728]
    refcols = match_to_list(lamp_lines['x'], rough_cols, plt=plt)
    # plt.show()

    def set_line_identity(lamp_lines, linelist, x, wav):
        # find the closest matching lamp_line
        ilamp = np.argmin(np.abs(lamp_lines['x'] - x))
        
        # find the closest matching row in the linelist
        ilist = np.argmin(np.abs(linelist['wav'] - wav))
        
        # reset values in lamp_lines
        lamp_lines[ilamp]['id'] = linelist[ilist]['id']
        lamp_lines[ilamp]['wav'] = linelist[ilist]['wav']
        # record ids and wavelengths of matched lines
    for col, wav in zip(refcols, refwavs):
        set_line_identity(lamp_lines, linelist, col, wav)

    # this routine finds the polynomial coefficients needed to transform between 
    # column numbers and wavelengths (and vice versa). Outlier rejection is included.
    def get_wavelength_solution(lamp_lines, order=4):
        wfit = np.isfinite(lamp_lines['wav'])
        
        # define the reverse mapping (wavelength to column)
        revcoeff = np.polyfit(lamp_lines['wav'][wfit], lamp_lines['x'][wfit], order)
    
        # define the forward mapping (column to wavelength)
        coeff = np.polyfit(lamp_lines['x'][wfit], lamp_lines['wav'][wfit], order)
        
        # check the fit for outliers
        fit_wav = np.polyval(coeff, lamp_lines['x'])
        wavres = fit_wav - lamp_lines['wav']
        lamp_lines['wavres'] = wavres
        sample = wavres[wfit]
        sample.sort()
        sample = sample[int(0.1 * sample.size) : int(0.9 * sample.size + 0.5)]    
        std = np.std(sample, ddof=1)
        w = wfit 
        w[wfit] = (np.abs(lamp_lines['wavres'][wfit]) < (5 * std))
        if w.sum() != lamp_lines.size:
            # re-fit with outliers rejected
            coeff, revcoeff = get_wavelength_solution(lamp_lines[w], order=order)
            
            # reset wavelength residuals using new coefficients
            fit_wav = np.polyval(coeff, lamp_lines['x'])
            wavres = fit_wav - lamp_lines['wav']
            lamp_lines['wavres'] = wavres
            
        lamp_lines['used'] = w
        return coeff, revcoeff
            
    def check_wavelength_solution(lamp_spec, lamp_lines, coeff):    
        wavs = col_to_wav(coeff, np.arange(lamp_spec.size))
        plt.plot(wavs, lamp_spec, c='g', lw=2)
        mark_peaks(plt, lamp_lines, 'wav')
        plt.colorbar(label='Residual ($\AA$)')
        plt.xlabel('Wavelength ($\AA$)')
        plt.ylabel('Counts')
        plt.xlim(min(wavs),max(wavs))
        plt.grid()
        
    def col_to_wav(coeff, cols):
        return np.polyval(coeff, cols)
    
    def wav_to_col(revcoeff, wavs):
        return np.polyval(revcoeff, wavs)
    
    # def mark_matched(lamp_lines):
    #     for line in lamp_lines:
    #         plt.axvline(line['wav'], ls='dotted', c='k')
    #         plt.text(line['wav'],max(lamp_ref['counts'])*.5,str(np.round(line['wav'],2)))
            # plt.show()
    # estimate a linear relation between column number and wavlength

    def mark_matched(refwavs):
        for line in refwavs:
            plt.axvline(line, ls='dotted', c='k',alpha=.5)
            t = plt.text(line,max(lamp_ref['counts'])*.9,str(np.round(line,2)),
                rotation=90,ha='center', va='center')
            t.set_bbox(dict(facecolor='white', edgecolor='white',alpha=.7))
            # plt.show()


    # def mark_matched(lamp_lines):
    #     # for line in lamp_lines:
    #     #     plt.axvline(line['wav'], ls='dotted', c='k')
    #     lamp_lines_n = len(lamp_lines)
    #     mark = range(0,lamp_lines_n,20)
    #     print (list(mark),mark)
    #     input()
    #     # print (mark)
    #     for i in mark:
    #         # print (i)
    #         # input()
    #         plt.axvline(lamp_lines['wav'][i], ls='dotted', c='k')
    #         plt.text(lamp_lines['wav'][i],max(lamp_ref['counts'])*.5,str(np.round(lamp_lines['wav'][i],2)))



    coeff, revcoeff = get_wavelength_solution(lamp_lines, order=order)
    
    # apply the wavelength solution to the column numbers
    wavs = col_to_wav(coeff, np.arange(lamp_spec.size))

    # print (wavs)
    # print (coeff)
    # print (lamp_spec)
    # input('Enter')
    

    # if display is True:
    #     # plot the initial wavelength calibrated spectrum
    #     plt.close()
    #     plt.figure(figsize=(15, 7))
    #     plt.plot(wavs, lamp_spec, c='g', lw=2, label='Arc spectrum')
    #     # plot the reference spectrum in red
    #     plt.plot(lamp_ref['wav'], lamp_ref['counts'], label='Reference', c='r')
    #     plt.xlim(wav1, wav2)
    #     plt.xlabel('Wavelength ($\AA$)'); plt.ylabel('Counts')
    #     mark_matched(lamp_lines)
    #     plt.legend(loc='best')
    #     plt.title('Polynomial fitting order: '+ str(order))
    #     plt.savefig(target_name+'_'+file_name+'_traced_arc.png')
    #     plt.show()

    # if display is False:
    #     # plot the initial wavelength calibrated spectrum
    #     plt.close()
    #     plt.figure(figsize=(15, 7))
    #     plt.plot(wavs, lamp_spec, c='g', lw=2, label='Arc spectrum')
    #     # plot the reference spectrum in red
    #     plt.plot(lamp_ref['wav'], lamp_ref['counts'], label='Reference', c='r')
    #     plt.xlim(wav1, wav2)
    #     plt.xlabel('Wavelength ($\AA$)'); plt.ylabel('Counts')
    #     mark_matched(lamp_lines)
    #     plt.legend(loc='best')
    #     plt.title('Polynomial fitting order: '+ str(order))
    #     plt.savefig(target_name+'_'+file_name+'_traced_arc.png')
    #     plt.close()


    # check for more matches in the range already fit
    def match_more(lamp_lines, linelist, order=4, tol=2):
        coeff, revcoeff = get_wavelength_solution(lamp_lines, order=order)
        wfit = np.isfinite(lamp_lines['wav'])
        minwav = lamp_lines['wav'][wfit].min()
        maxwav = lamp_lines['wav'][wfit].max()
        
        xmin = lamp_lines['x'][wfit].min()
        xmax = lamp_lines['x'][wfit].max()
        
        w = (lamp_lines['x'] > xmin) & (lamp_lines['x'] < xmax)
        for line in lamp_lines[w]:
            rough_wav = col_to_wav(coeff, line['x'])
            refwav = match_to_list(linelist['wav'], [rough_wav], tol=tol)
            if np.isfinite(refwav):
                #print(f'matched column {line["x"]:.1f} to wavelength {refwav[0]}')
                set_line_identity(lamp_lines, linelist, line['x'], refwav)
    

    # match_more(lamp_lines, linelist, order=order)
    # # re-fit with a higher order

    # coeff, revcoeff = get_wavelength_solution(lamp_lines, order=4)
    # wavs = col_to_wav(coeff, np.arange(lamp_spec.size))

    # # check the supposed matches
    # col1, col2 = np.polyval(revcoeff, [wav1, wav2])
    # refcols = match_to_list(lamp_lines['x'], rough_cols, plt=None)
    # for col, wav in zip(refcols, refwavs):
    #     set_line_identity(lamp_lines, linelist, col, wav)
        
    # # auto-match more lines
    # match_more(lamp_lines, linelist)
    
    # # re-fit
    # coeff, revcoeff = get_wavelength_solution(lamp_lines, order=4)
    # check_wavelength_solution(lamp_spec, lamp_lines, coeff)

    match_more(lamp_lines, linelist, order=order)
    # re-fit with a higher order

    coeff_more, revcoeff_more = get_wavelength_solution(lamp_lines, order=4)
    wavs_more = col_to_wav(coeff_more, np.arange(lamp_spec.size))

    # check the supposed matches
    col1, col2 = np.polyval(revcoeff_more, [wav1, wav2])
    refcols = match_to_list(lamp_lines['x'], rough_cols, plt=None)
    for col, wav in zip(refcols, refwavs):
        set_line_identity(lamp_lines, linelist, col, wav)
        
    # auto-match more lines
    match_more(lamp_lines, linelist)
    
    # re-fit
    coeff_more, revcoeff_more = get_wavelength_solution(lamp_lines, order=4)
    check_wavelength_solution(lamp_spec, lamp_lines, coeff_more)


    if display is True:
        # plt.xlim(wave_min,wave_max)
        plt.savefig(target_name+'_'+file_name+'_wavelength_sol.png')
        plt.show()

    if display is False:
        # plt.xlim(wave_min,wave_max)
        plt.savefig(target_name+'_'+file_name+'_wavelength_sol.png')
        plt.close()

    w = lamp_lines['used']
    std = np.std(lamp_lines['wavres'][w], ddof=1)
    print(f'STD of wavelength residual is {std:0.2} Angstrom')
    
    # f = open('ref_cuar_gr7_162.dat','w')
    # for i in range(len(wavs)):
    #     f.write('%f %f\n' % (wavs[i],lamp_spec[i]))
    # f.close()


    # Plot residuals
    # plt.scatter(lamp_lines['wav'][w], lamp_lines['wavres'][w]) #, c=np.abs(lamp_lines['wavres'][w]))
    # plt.scatter(lamp_lines['wav'][~w], lamp_lines['wavres'][~w], marker='x') #c=np.abs(lamp_lines['wavres'][~w]))
    # plt.axhspan(-std, std, color='k', alpha=0.1)
    # plt.xlabel('Wavelength ($\AA$)'); plt.ylabel('Counts'); plt.grid()
    # plt.show()

    pixel_scale = (max(wavs) - min(wavs))/len(lamp_spec)

    if display is True:
        # plot the initial wavelength calibrated spectrum
        plt.close()
        plt.figure(figsize=(15, 7))
        plt.plot(wavs, lamp_spec, c='g', lw=2, label='Arc spectrum')
        # plot the reference spectrum in red
        plt.plot(lamp_ref['wav'], lamp_ref['counts'], label='Reference', c='r')
        plt.xlim(wav1, wav2)
        plt.xlabel('Wavelength ($\AA$)'); plt.ylabel('Counts')
        mark_matched(refwavs)
        plt.legend(loc='best')
        plt.title('Polynomial fitting order: '+ str(order)+ ', \
            STD of wavelength residual: '+str(np.round(std,2)) +' $\AA$,'\
            + ' Pixel scale = ' + str(np.round(pixel_scale,2)) +' $\AA$/pix')
        plt.savefig(target_name+'_'+file_name+'_traced_arc.png')
        plt.show()

    if display is False:
        # plot the initial wavelength calibrated spectrum
        plt.close()
        plt.figure(figsize=(15, 7))
        plt.plot(wavs, lamp_spec, c='g', lw=2, label='Arc spectrum')
        # plot the reference spectrum in red
        plt.plot(lamp_ref['wav'], lamp_ref['counts'], label='Reference', c='r')
        plt.xlim(wav1, wav2)
        plt.xlabel('Wavelength ($\AA$)'); plt.ylabel('Counts')
        mark_matched(refwavs)
        plt.legend(loc='best')
        plt.title('Polynomial fitting order: '+ str(order)+ \
            ', STD of wavelength residual: '+str(np.round(std,2)) +' $\AA$,'\
            + ' Pixel scale = ' + str(np.round(pixel_scale,2)) +' $\AA$/pix')
        plt.savefig(target_name+'_'+file_name+'_traced_arc.png')
        plt.close()


    #

    spec_smooth_sum = smooth(spec_sum,3)
    if onedspec_optimal is not None:
        spec_smooth_optimal = smooth(spec_optimal,3)



    if display is True:
        plt.plot(wavs_more, spec_smooth_sum,linewidth=1.,label='Summed trace')
        if onedspec_optimal is not None:
            plt.plot(wavs_more, spec_smooth_optimal,linewidth=1.,label='Optimal trace',color='orange')
        plt.title('Reduced spectrum of '+spec_file_name+': ' +target_name)
        plt.xlabel('Wavelength $\AA$')
        plt.ylabel('Counts')
        plt.legend(loc='best')
        plt.savefig(target_name+'_'+spec_file_name+'_reduced_spectrum.png')
        plt.show()

    if display is False:
        plt.plot(wavs_more, spec_smooth_sum,linewidth=1.,label='Summed trace')
        if onedspec_optimal is not None:
            plt.plot(wavs_more, spec_smooth_optimal,linewidth=1.,label='Optimal trace',color='orange')
        plt.title('Reduced spectrum of '+spec_file_name+': ' +target_name)
        plt.xlabel('Wavelength $\AA$')
        plt.ylabel('Counts')
        plt.legend(loc='best')
        plt.savefig(target_name+'_'+spec_file_name+'_reduced_spectrum.png')
        plt.close()

    # fname = k.split('.')[0]
    f = open(target_name+'_'+spec_file_name+'_reduced.dat','w')
    if onedspec_optimal is not None:
        f.write('wavelength spectrum_sum 3smoothed_sum spectrum_optimal 3smoothed_optimal\n')
        for ii in range(len(wavs)):
            f.write('%g %g %g %g %g\n' % (wavs_more[ii], spec_sum[ii],spec_smooth_sum[ii],spec_optimal[ii],spec_smooth_optimal[ii]))
        f.close()
    else:
        f.write('wavelength spectrum_sum 3smoothed_sum\n')
        for ii in range(len(wavs)):
            f.write('%g %g %g\n' % (wavs_more[ii], spec_sum[ii],spec_smooth_sum[ii]))
        f.close()

    return interact

    # if display is True:
    #     plt.plot(wavs, spec_smooth_sum,linewidth=1.,label='Summed trace')
    #     if onedspec_optimal is not None:
    #         plt.plot(wavs, spec_smooth_optimal,linewidth=1.,label='Optimal trace',color='orange')
    #     plt.title('Reduced spectrum of '+spec_file_name+': ' +target_name)
    #     plt.xlabel('Wavelength $\AA$')
    #     plt.ylabel('Counts')
    #     plt.legend(loc='best')
    #     plt.savefig(target_name+'_'+spec_file_name+'_reduced_spectrum.png')
    #     plt.show()

    # if display is False:
    #     plt.plot(wavs, spec_smooth_sum,linewidth=1.,label='Summed trace')
    #     if onedspec_optimal is not None:
    #         plt.plot(wavs, spec_smooth_optimal,linewidth=1.,label='Optimal trace',color='orange')
    #     plt.title('Reduced spectrum of '+spec_file_name+': ' +target_name)
    #     plt.xlabel('Wavelength $\AA$')
    #     plt.ylabel('Counts')
    #     plt.legend(loc='best')
    #     plt.savefig(target_name+'_'+spec_file_name+'_reduced_spectrum.png')
    #     plt.close()

    # # fname = k.split('.')[0]
    # f = open(target_name+'_'+spec_file_name+'_reduced.dat','w')
    # if onedspec_optimal is not None:
    #     f.write('wavelength spectrum_sum 3smoothed_sum spectrum_optimal 3smoothed_optimal\n')
    #     for ii in range(len(wavs)):
    #         f.write('%g %g %g %g %g\n' % (wavs[ii], spec_sum[ii],spec_smooth_sum[ii],spec_optimal[ii],spec_smooth_optimal[ii]))
    #     f.close()
    # else:
    #     f.write('wavelength spectrum_sum 3smoothed_sum\n')
    #     for ii in range(len(wavs)):
    #         f.write('%g %g %g\n' % (wavs[ii], spec_sum[ii],spec_smooth_sum[ii]))
    #     f.close()

    # return interact


def flux_callibration(standard_reduced_spec,standard_name,science_spec,display):
    # https://www.eso.org/sci/observing/tools/standards/spectra/
    file_name = science_spec.split('.')[0].split('_reduced')[0]
    file_name_std = standard_reduced_spec.split('.')[0].split('_reduced')[0]

    standard_spec = np.loadtxt(standard_reduced_spec,skiprows=1,unpack=True)
    
    if len(standard_spec) == 3:
        w_stand, sum_stand = standard_spec[0], standard_spec[1]

    if len(standard_spec) == 5:
        w_stand, sum_stand, opt_stand = standard_spec[0], standard_spec[1], standard_spec[3]
    # w, f = np.loadtxt(standard_reduced_spec,skiprows=1,usecols=(0,1),unpack=True)

    data =  np.loadtxt(science_spec,skiprows=1,unpack=True)
    if len(data) == 3:
        optimal = None
        waves_cor, sum_spec = data[0], data[1]

    if len(data) == 5:
        optimal = True
        waves_cor, sum_spec, opt_spec = data[0], data[1], data[3]

    # Download the data from ESO
    dtype = []
    dtype.append( ('wav', float) )
    dtype.append( ('flux', float) ) # units are ergs/cm/cm/s/A * 10**16
    dtype.append( ('eflux', float) )
    dtype.append( ('dlam', float) )
    # calspec = np.genfromtxt('ftp://ftp.eso.org/pub/stecf/standards/okestan/ffeige67.dat', dtype=dtype)
    # calspec = np.genfromtxt('fltt3218.dat', dtype=dtype)

    def response(w,f):
        if standard_name == 'ltt3218':
            calspec = np.genfromtxt('ftp://ftp.eso.org/pub/stecf/standards/ctiostan/fltt3218.dat',dtype=dtype)
        if standard_name == 'eg21':
            calspec = np.genfromtxt('ftp://ftp.eso.org/pub/stecf/standards/ctiostan/feg21.dat',dtype=dtype)
        if standard_name == 'ltt377':   #CD-34d241 was mistakenly named LTT 377 in Hamuy et al (1992 & 1994)
            calspec = np.genfromtxt('ftp://ftp.eso.org/pub/stecf/standards/ctiostan/fcd_34d241.dat',dtype=dtype)
        if standard_name == 'feige110':  
            calspec = np.genfromtxt('ftp://ftp.eso.org/pub/stecf/standards/ctiostan/ffeige110.dat',dtype=dtype)
        if standard_name == 'cd-32-9927':  
            calspec = np.genfromtxt('ftp://ftp.eso.org/pub/stecf/standards/ctiostan/fcd32d9927.dat',dtype=dtype)
        if standard_name == 'ltt7379':  
            calspec = np.genfromtxt('ftp://ftp.eso.org/pub/stecf/standards/ctiostan/fltt7379.dat',dtype=dtype)
        # fit a spline to the tabulated spectrum
        t = np.arange(calspec['wav'][1], calspec['wav'][-2], np.int(np.median(calspec['dlam'])))
        # print (calspec['dlam'])
        stdflux = interpolate.InterpolatedUnivariateSpline(calspec['wav'], calspec['flux'])
        std_spec = f
        waves = w
        ratios = std_spec / stdflux(waves)
        w = (waves> calspec['wav'].min()) \
            & (waves < calspec['wav'].max()) \
            & (np.abs(waves - 7650) > 70) \
            & (np.abs(waves - 6900) > 40) \
            & (np.abs(waves - 6563) > 40) \
        # fit a spline to the ratios to determine the response function
        t = waves[w][1:-2:25]
        respfn = LSQUnivariateSpline(waves[w], ratios[w], t)
        return respfn, waves, ratios, w, std_spec, calspec

    try:
        if len(data) == 3 :
                respfn, waves, ratios, w, std_spec, calspec = response(w_stand, sum_stand)
                respfn_optimal = None
        elif (len(data) == 5) and (len(standard_spec) == 5):
            respfn, waves, ratios, w, std_spec, calspec = response(w_stand, sum_stand)
            respfn_optimal, waves_optimal, ratios_optimal, w_optimal, std_spec_optimal, calspec_optimal = response(w_stand, opt_stand)

        f = open(science_spec.split('_reduced.dat')[0]+'_fluxcal.dat','w')
        if respfn_optimal is not None:
            f.write('wavelength spec_sum spec_optimal\n')
            for ii in range(7,len(waves_cor)-3):
                f.write('%g %g %g\n' % (waves_cor[ii], sum_spec[ii] / respfn(waves_cor[ii]),opt_spec[ii] / respfn_optimal(waves_cor[ii])))
            f.close()
        else:
            f.write('wavelength spec_sum\n')
            for ii in range(7,len(waves_cor)-3):
                f.write('%g %g\n' % (waves_cor[ii], sum_spec[ii] / respfn(waves_cor[ii])))
            f.close()
        
        #
        #Plot response function
        #
        if respfn_optimal is not None:
            plt.subplot(211)
        plt.title('Calibrating the standard:' + standard_reduced_spec.split('_reduced.dat')[0])
        plt.plot(waves[w], ratios[w], 'ro',ms=1.,label='Summed spectrum')
        xwav = np.linspace(waves[w][1], waves[w][-1], 1000)
        plt.plot(xwav, respfn(xwav))
        plt.xlabel('Wavelength ($\AA$)')
        plt.ylabel('Response Function $\\left(\\frac{Counts}{F_{\lambda}}\\right)$')
        plt.legend(loc=1)

        if respfn_optimal is not None:
            plt.subplot(212)
            # plt.title('Calibrated standard: ' + standard_reduced_spec.split('_reduced.dat')[0])
            plt.plot(waves_optimal[w_optimal], ratios_optimal[w_optimal], 'ro',ms=1.,label='Optimal spectrum')
            xwav = np.linspace(waves_optimal[w_optimal][1], waves_optimal[w_optimal][-1], 1000)
            plt.plot(xwav, respfn_optimal(xwav))
            plt.xlabel('Wavelength ($\AA$)')
            plt.ylabel('Response Function $\\left(\\frac{Counts}{F_{\lambda}}\\right)$')
            plt.legend(loc=1)
        
        plt.savefig(file_name_std+'_fluxcal_response.png')
        if (display is True):
            plt.show()
        plt.clf()

        #
        #Plot calibrated standard
        #
        if respfn_optimal is not None:
            plt.subplot(211)

        plt.title('Calibrated standard: ' + standard_reduced_spec)
        plt.plot(calspec['wav'][7:-3], calspec['flux'][7:-3], label='Tabulated (published) Spectrum',color='orange')
        plt.plot(waves[7:-3], std_spec[7:-3] / respfn(waves[7:-3]), label='Summed extracted Spectrum')
        plt.xlabel('Wavelength ($\AA$)')
        plt.ylabel('$F_{\lambda} \: \left(\\frac{ergs}{cm^2\, s\, \AA} \: \\times 10^{16}\\right)$')
        plt.yscale('log')
        plt.legend(loc=1)

        if respfn_optimal is not None:
            plt.subplot(212)
            # plt.title('Calibrated standard: ' + standard_reduced_spec)
            plt.plot(calspec_optimal['wav'][7:-3], calspec_optimal['flux'][7:-3], label='Tabulated (published) Spectrum',color='orange')
            plt.plot(waves_optimal[7:-3], std_spec_optimal[7:-3] / respfn_optimal(waves_optimal[7:-3]), label='Optimal extracted spectrum')
            plt.xlabel('Wavelength ($\AA$)')
            plt.ylabel('$F_{\lambda} \: \left(\\frac{ergs}{cm^2\, s\, \AA} \: \\times 10^{16}\\right)$')
            plt.yscale('log')
            plt.legend(loc=1)

        plt.savefig(file_name_std+'_fluxcal_standard.png')
        if display is True:
            plt.show()
        plt.clf()
        #
        # Plot flux calibrated spectra
        #
        plt.figure(figsize=(10, 5))
        plt.title('Flux calibrated: '+science_spec.split('_reduced.dat')[0])
        if optimal is not None:
            plt.plot(waves_cor[7:-3], opt_spec[7:-3] / respfn_optimal(waves_cor[7:-3]),label='Optimal trace',color='orange')
        plt.plot(waves_cor[7:-3], sum_spec[7:-3] / respfn(waves_cor[7:-3]),label='Summed trace')
        plt.legend(loc='best')
        plt.xlabel('Wavelength ($\AA$)')
        plt.ylabel('$F_{\lambda} \: \left(\\frac{ergs}{cm^2\, s\, \AA} \: \\times 10^{16}\\right)$')
        plt.savefig(science_spec.split('_reduced.dat')[0]+'_fluxcal.png')
        if display is True:
            plt.show()
 
    except NameError:
        print ('Standard not found, add standard eso url to script')
        print ('File: '+science_spec)



def combine_dat(data_files,out_file,spectrum_keyword,wavelength_keyword,combine_method='median',display=True,sigma=3,accuray=2):
    '''
    data_files: list / str
        list of files to be combined
    spectrum_keyword: str
        Which spectrum to combine from the reduced data files or fluxcal data files
    wavelength_keyword: str
        Identifies wavlength keyword from reduced data files or fluxcal data files
    combine_method: str
        Method to combine spectra fluxes.  Default median, options are:
        - median
        - mean
        - sum
    sigma:  int
        sigma_clipping, default sigma =  3
    accuray:  int
        depending on resolution, may need to increase for higher resolution --> integer
        default = 2
    output:
    Save file as space delimited with median wavelength and combined spectrum)
    '''
    n = len(data_files)
    frames = ['df'+str(i) for i in range(n)]
    d = [pd.read_csv(data_files[i],delimiter=' ') for i in range(n)]
    #
    pix_scale = np.round(max([(i[wavelength_keyword].shift(1)-i[wavelength_keyword]).abs().mean() for i in d]),accuray)

    lam_min_id, lam_min = np.array([i[wavelength_keyword].min() for i in d]).argmax(), max([i[wavelength_keyword].min() for i in d])
    lam_max_id, lam_max = np.array([i[wavelength_keyword].max() for i in d]).argmin(), min([i[wavelength_keyword].max() for i in d])

    # print (lam_min_id, lam_min,lam_max_id, lam_max)


    df3 = pd.DataFrame()
    names_fluxes = np.empty(n*2,dtype=object)
    names_waves = np.empty(n*2,dtype=object)
    names_fluxes_clipped = []
    names_fluxes_clipped = []

    for i in range(n):
        d_lam_min = (d[i][wavelength_keyword] >= lam_min - 1*pix_scale) & (d[i][wavelength_keyword] <= lam_min + 1*pix_scale)
        min_index = abs(d[i][d_lam_min][wavelength_keyword] - lam_min).idxmin()
        d_lam_max = (d[i][wavelength_keyword] >= lam_max - 1*pix_scale) & (d[i][wavelength_keyword] <= lam_max + 1*pix_scale)
        max_index = abs(d[i][d_lam_max][wavelength_keyword] - lam_max).idxmin()

        name_flux = 'flux '+ str(i)
        name_wav = 'wave ' + str(i)
    
        sig = sigma_clip(d[i][spectrum_keyword][min_index:max_index],\
            sigma = sigma,\
            cenfunc=np.mean,\
            masked=True)
        mask_sig = d[i][spectrum_keyword][min_index:max_index] == sig
        wave_match_clipped = d[i][min_index:max_index][mask_sig.fillna(False)]
        raw_flux =  d[i][spectrum_keyword][min_index:max_index]
        raw_wave =  d[i][wavelength_keyword][min_index:max_index]
        
        if i == 0:
            names_fluxes[i] = name_flux
            names_waves[i] = name_wav
            names_fluxes[i+1] = name_flux+'-clipped'
            names_waves[i+1] = name_wav+'-clipped'

        else:
            names_fluxes[2*i] = name_flux
            names_waves[2*i] = name_wav
            names_fluxes[2*i+1] = name_flux+'-clipped'
            names_waves[2*i+1] = name_wav+'-clipped'

        
        name_wav_clipped = name_wav+'-clipped'
        name_flux_clipped = name_flux+'-clipped'
        name_wav = pd.DataFrame({name_wav: raw_wave})
        name_flux = pd.DataFrame({name_flux: raw_flux})
        name_wav_clipped = pd.DataFrame({name_wav_clipped: wave_match_clipped[wavelength_keyword]})
        name_flux_clipped = pd.DataFrame({name_flux_clipped: wave_match_clipped[spectrum_keyword]})
        
        if i == 0:
            df3 = pd.concat([name_wav,name_flux,name_wav_clipped,name_flux_clipped],axis=1)
        else:
            df3 = pd.concat([df3,name_wav,name_flux,name_wav_clipped,name_flux_clipped],axis=1)

    df3 = df3.apply(lambda x: pd.Series(x.dropna().to_numpy())) 


    clipped_fluxes = list(filter(lambda x: ('clipped' in x) and ('flux' in x),df3.columns))
    clipped_waves = list(filter(lambda x: ('clipped' in x) and ('wave' in x),df3.columns))

    zipped_fluxes = list(zip(*[df3[c] for c in clipped_fluxes]))
    zipped_waves = list(zip(*[df3[c] for c in clipped_waves]))

    df3['median waves'] = [np.median(i) for i in zipped_waves]
    col_name = 'clipped '+combine_method+' fluxes'
    if combine_method == 'mean':
        df3[col_name] = [np.mean(i) for i in zipped_fluxes]
    elif combine_method == 'median':
        df3[col_name] = [np.median(i) for i in zipped_fluxes]
    elif combine_method == 'sum':
        df3[col_name] = [np.sum(i) for i in zipped_fluxes]
    

    df3[['median waves',col_name]].to_csv(out_file+'_'+combine_method+'.dat',sep=' ',na_rep=np.nan, index=False, header=True)


    if display is True:

        fluxes = list(filter(lambda x: ('clipped' not in x) and ('flux' in x),df3.columns))
        waves = list(filter(lambda x: ('median' not in x) and ('clipped' not in x) and ('wave' in x),df3.columns))
        [plt.plot(df3[list(clipped_waves)[i]],df3[list(clipped_fluxes)[i]],\
            label='{}, std = {}'.format(list(clipped_fluxes)[i],df3[list(clipped_fluxes)[i]].std().\
                round(2))) for i in range(n)] 
        plt.plot(df3['median waves'],df3[col_name],label='{} combined, std = {}'.format(combine_method,df3[col_name].std().round(2)))
        plt.legend(loc=1)
        plt.xlabel('Median Wavelength ($\AA$)')
        plt.savefig(out_file+'_'+combine_method+'.png')
        plt.show()

    if display is False:

        fluxes = list(filter(lambda x: ('clipped' not in x) and ('flux' in x),df3.columns))
        waves = list(filter(lambda x: ('median' not in x) and ('clipped' not in x) and ('wave' in x),df3.columns))
        [plt.plot(df3[list(clipped_waves)[i]],df3[list(clipped_fluxes)[i]],\
            label='{}, std = {}'.format(list(clipped_fluxes)[i],df3[list(clipped_fluxes)[i]].std().\
                round(2))) for i in range(n)] 
        plt.plot(df3['median waves'],df3[col_name],label='{} combined, std = {}'.format(combine_method,df3[col_name].std().round(2)))
        plt.legend(loc=1)
        # # # # plt.gca().set_ylim(bottom=0)
        plt.xlabel('Median Wavelength ($\AA$)')
        plt.savefig(out_file+'.png')    


def normalise(file,wavelength_keyword,spectrum_keyword,exclude_range,interactive,display,out_file,poly_order=8):
    '''
    file: str
        Data file (ascii) that contains the spectrum you would like to normalise
    wavelength_keyword: str or int
        str:  When your data file have column names, identify the correct column name for
              the wavelength 
        int:  When your data file has no column name identify the wavelength column number,
              start counting form 0
    spectrum_keyword: str or int
        Same as wavelength_keyword but for flux (i.e the spectrum)
    display: True / False
        True:   Plot will be displayed and saved as png
        False:  Plot will not be displayed, but saved as png
    out_file: str
        The name of the output files, but WITHOUT EXTENSION.  Two files will be generated.
        A plot as well as a data file with normalised spectrum.
    order:  int
        Polynomial that will be fit to spectrum, default 8
    '''
    df = pd.read_csv(file,delimiter=' ')
    if (type(wavelength_keyword) is not str) or (type(spectrum_keyword) is not str):
        w = df.columns[wavelength_keyword]
        f = df.columns[spectrum_keyword]
    else:
        w = wavelength_keyword
        f =  spectrum_keyword
    pos = df[f] > 0
    df = df[pos]

    df_out = pd.DataFrame()

    c = np.polyfit(x=df[w],y=df[f], deg=poly_order)
    p = np.poly1d(c)
    # fit = np.polyval(c,df[w])
    df['norm'] = df[f]/p(df[w])

    def normal_plot(df,title):
        fig = plt.figure(figsize=(7,6))
        plt.subplot(211)
        plt.title(title)
        plt.plot(df[w],df[f],label='Spectrum',color='royalblue',alpha=.8)
        plt.plot(df[w],p(df[w]),label='Fit: poly order = {}'.format(poly_order),color='k',ls='--')
        plt.xlabel('Wavelength ($\AA$)')
        plt.legend(loc=1)
        plt.subplot(212)
        plt.plot(df[w],df['norm'],label='Normalised spectrum',color='royalblue',alpha=.8)
        plt.xlabel('Wavelength ($\AA$)')
        plt.legend(loc=1)
        plt.tight_layout()

    def normal_plot_refit(df,p,df_new,p_new,title):
        fig = plt.figure(figsize=(7,6))
        plt.subplot(211)
        plt.title(title)
        plt.plot(df[w],df[f],label='Excluded spectrum',color='royalblue',alpha=.2)
        plt.plot(df_new[w],df_new[f],label='Spectrum',color='royalblue')
        plt.plot(df[w],p(df[w]),label='Fit: poly order = {}'.format(poly_order),color='k',alpha=.8)
        plt.plot(df[w],p_new(df[w]),label='Adjusted poly fit',color='red',ls='--')
        plt.xlabel('Wavelength ($\AA$)')
        plt.legend(loc=1)
        plt.subplot(212)
        plt.plot(df[w],df['norm'],label='Normalised spectrum',color='k')
        plt.plot(df[w],df[f]/p_new(df[w]),label='Adjusted normalised spectrum',color='royalblue',alpha=.8)
        plt.xlabel('Wavelength ($\AA$)')
        plt.legend(loc=1)
        plt.tight_layout()


    k = 0

    if (interactive is False) and (exclude_range is None):
        df_out[w] = df[w]
        df_out['norm'] = df['norm']
        df_out.to_csv(out_file+'_norm.dat',sep=' ',na_rep=np.nan, index=False, header=True)
        normal_plot(df=df,title='Normalisation of {}'.format(file))

        if display is False:
            plt.savefig(out_file+'_norm.png')
        if display is True:
            plt.savefig(out_file+'_norm.png')
            plt.show()

    if (interactive is False) and (exclude_range is not None):
        mask_ranges = [(df[w] >= exclude_range[i]) & (df[w] <= exclude_range[i+1]) for i in range(0,len(exclude_range)-1,2)]
        df_ranges = df
        for i in mask_ranges:
            df_ranges = df_ranges[~i]
        
        c_new = np.polyfit(x=df_ranges[w],y=df_ranges[f], deg=poly_order)
        p_new = np.poly1d(c_new)

        normal_plot_refit(df,p,df_ranges,p_new,title='Normalisation of {}'.format(file))

        df_out[w] = df[w]
        df_out['norm'] = df['norm']
        df_out.to_csv(out_file+'_norm.dat',sep=' ',na_rep=np.nan, index=False, header=True)

        if display is False:
            plt.savefig(out_file+'_norm.png')
        if display is True:
            plt.savefig(out_file+'_norm.png')
            plt.show()



    if interactive is True:
        k = 0
        def interactive_ranges(mask_ranges,exclude_range,k):
            if (mask_ranges is None) or (exclude_range is not None):
                if (mask_ranges is None) and (exclude_range is None):
                    normal_plot(df=df,title='Which ranges to exclude?')
                    plt.show()
                    ranges = np.array([float(i) for i in np.array(input('Insert ranges (e.g. 4000 4050 5080 6000): ').split(' '))])
                    mask_ranges = [(df[w] >= ranges[i]) & (df[w] <= ranges[i+1]) for i in range(0,len(ranges)-1,2)]
                elif exclude_range is not None:
                    mask_ranges = [(df[w] >= exclude_range[i]) & (df[w] <= exclude_range[i+1]) for i in range(0,len(exclude_range)-1,2)]
                    df_ranges = df
                    for i in mask_ranges:
                        df_ranges = df_ranges[~i]
                    
                    c_new = np.polyfit(x=df_ranges[w],y=df_ranges[f], deg=poly_order)
                    p_new = np.poly1d(c_new)
                    
                    normal_plot_refit(df,p,df_ranges,p_new,title='Excluded spectrum, exclude_range more regions?')
                    plt.show()
                    cont = input('Exclude another region (y/n)? ')
                    k += 1
                    exclude_range = None
                    if cont == 'n':
                        plt.savefig(out_file+'_norm.png')
                        df_out[w] = df[w]
                        df_out['norm'] = df[f]/p_new(df[w])
                        df_out.to_csv(out_file+'_norm.dat',sep=' ',na_rep=np.nan, index=False, header=True)
                        sys.exit()

            if (mask_ranges is not None):
                if k > 0:
                    ranges = np.array([float(i) for i in np.array(input('Insert ranges (e.g. 4000 4050 5080 6000): ').split(' '))])
                    mask_ranges_old = mask_ranges
                    mask_ranges = [(df[w] >= ranges[i]) & (df[w] <= ranges[i+1]) for i in range(0,len(ranges)-1,2)]
                    mask_ranges.extend(mask_ranges_old)
                else:
                    mask_ranges_old = mask_ranges
                    mask_ranges = [(df[w] >= ranges[i]) & (df[w] <= ranges[i+1]) for i in range(0,len(ranges)-1,2)]
                    mask_ranges.extend(mask_ranges_old)

                
            df_ranges = df
            for i in mask_ranges:
                df_ranges = df_ranges[~i]
            
            c_new = np.polyfit(x=df_ranges[w],y=df_ranges[f], deg=poly_order)
            p_new = np.poly1d(c_new)
            
            normal_plot_refit(df,p,df_ranges,p_new,title='Excluded spectrum')
            plt.show()
            cont = input('Exclude another Region (y/n)? ')
            if cont == 'y':
                interactive = True
                exclude_range = None
                k += 1
                return interactive, mask_ranges, k, exclude_range
            else:
                normal_plot_refit(df,p,df_ranges,p_new,title='Excluded spectrum') 
                plt.savefig(out_file+'_norm.png')
                interactive = False
                df_out[w] = df[w]
                df_out['norm'] = df[f]/p_new(df[w])
                df_out.to_csv(out_file+'_norm.dat',sep=' ',na_rep=np.nan, index=False, header=True)
            return interactive, mask_ranges, k, exclude_range

        interactive, mask_ranges, k, exclude_range = interactive_ranges(mask_ranges=None,exclude_range=exclude_range,k=k)   

        while interactive is True:
            interactive, mask_ranges, k, exclude_range = interactive_ranges(mask_ranges=mask_ranges,exclude_range=exclude_range,k=k)

def smoothing(file,wavelength_keyword,spectrum_keyword,out_file,window=3,sigma=None,display=False):
    '''
    file: str
        Data file (ascii) that contains the spectrum you would like to normalise
    wavelength_keyword: str or int
        str:  When your data file have column names, identify the correct column name for
              the wavelength 
        int:  When your data file has no column name identify the wavelength column number,
              start counting form 0
    spectrum_keyword: str or int
        Same as wavelength_keyword but for flux (i.e the spectrum)
    display: True / False
        True:   Plot will be displayed and saved as png
        False:  Plot will not be displayed, but saved as png
    out_file: str
        The name of the output files, but WITHOUT EXTENSION.  Two files will be generated.
        A plot as well as a data file with normalised spectrum.
    window: int
        Window for smoothing, default = 3
    sigma: None / int
        If sigma is not None (default), spectrum will be sigma clipped withing window range
    '''

    df = pd.read_csv(file,delimiter=' ')
    if (type(wavelength_keyword) is not str) or (type(spectrum_keyword) is not str):
        w = df.columns[wavelength_keyword]
        f = df.columns[spectrum_keyword]
    else:
        w = wavelength_keyword
        f =  spectrum_keyword
    pos = df[f] > 0
    df = df[pos]

    smoothed = df.rolling(window=window).mean()
    smoothed_std = df.rolling(window=window).std()

    df_out = pd.DataFrame()

    if sigma is None:
        plt.plot(df[w],df[f],label='Spectrum')
        plt.plot(smoothed[w],smoothed[f],label='Smoothed - window = {}'.format(window))
        plt.legend(loc=1)
        plt.xlabel('Wavelength ($\AA$)')
        plt.title('Smoothed spectrum of {}'.format(file))

        df_out = smoothed
        df_out.to_csv(out_file+'_smoothed.dat',sep=' ',na_rep=np.nan, index=False, header=True)

        if display is True:
            plt.savefig(out_file+'_smoothed.png')
            plt.show()
        else: 
            plt.savefig(out_file+'_smoothed.png')


    if sigma is not None:

        sig = sigma_clip(smoothed_std,\
            sigma = sigma,\
            cenfunc = np.median,\
            masked=True)
        mask_sig = smoothed_std == sig
        df_clipped = smoothed[mask_sig.fillna(False)]

        df_out = df_clipped
        df_out.to_csv(out_file+'_smoothed_clipped.dat',sep=' ',na_rep=np.nan, index=False, header=True)

        fig = plt.figure(figsize=(8,6))

        plt.subplot(211)
        plt.title('Sigmal clipped, smoothed spectrum of {}'.format(file))
        plt.plot(df[w],df[f],label='Spectrum')
        plt.plot(smoothed[w],smoothed[f],label='Smoothed - window = {}'.format(window))
        plt.legend(loc=1)

        plt.subplot(212)
        plt.plot(df_clipped[w],df_clipped[f],label='Smoothed {}$\sigma$ clipped'.format(sigma))
        plt.legend(loc=1)
        plt.xlabel('Wavelength ($\AA$)')
        if display is True:
            plt.savefig(out_file+'_smoothed_clipped.png')
            plt.show()
        else: 
            plt.savefig(out_file+'_smoothed_clipped.png')


def fwhm_fits(file,row1,row2,column1,column2):       # pixel space
    cwd = os.getcwd()
    data_dir = cwd
    image = fits.getdata(os.path.join(data_dir,file))
    spec = np.sum(image[row1:row2, column1: column2], axis=0)  # 6965 A - closest to Halpha
    peaks, _ = find_peaks(spec, height=max(spec))

    results_half = peak_widths(spec, peaks, rel_height=0.5)
    # print (results_half)
    # print ('FWHM = ',results_half[0])
    plt.plot(spec)
    plt.plot(peaks, spec[peaks], "x")
    plt.hlines(*results_half[1:], color="C2")
    plt.text(results_half[3],results_half[1],'FWHM: ' +str(np.round(results_half[0][0],2)))
    # plt.show()

# fwhm_fits(file='a2241016.fits',row1=64,row2=73,column1=750,column2=820)




