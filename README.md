# Titus_Saures_Rex_Spectroscopic_Routines
This package includes a number of Python routines for spectroscopic reductions for long-slit spectra

### 1.	Functionality:

The pipeline uses native Python routines for the following tasks:
   * Produce master Bias
   * Produce master Flat
   * Apply calibrations:
       * Cosmic ray removal
       * Flat fielding
       * Bias subtraction
   * Trim images
   * Trace spectrum
   * Extract spectrum:
       * Sky subtraction
       * Sum trace
       * Optimal extraction
   * Wavelength calibration
   * Flux calibration
   * Smoothing
   * Combining 1D spectra
   * Normalise 1D spectra

### 2.  Dependancies (Python 3):

  *	Numpy
  *	os
  *	Astropy
  *	Specutils (See [installation](https://specutils.readthedocs.io/en/latest/))
  *	Matplotlib
  *	Scipy
  *	Lacosmic (See [installation](https://lacosmic.readthedocs.io/en/latest/))

All the dependancies can easily be installed using either pip or conda.

### 3.  The pipeline consists of two scripts:
  *	spectroscopic_routines.py
  *	spec.py

#### <ins>spectroscopic_routines.py:</ins> 
This is the script containing all the functions for the reductions.  You should not need to edit the script, but it contains the following functions:
  *	image_header():  Extracts value / parameter for specified header keyword
  *	trim():  trims images
  *	master_bias(): produce master bias
  *	master_flat(): produce master flat
  *	appy_calibrations():  Applies flat fielding, bias subtraction, cosmic ray removal
  *	twodsky():  Produces a 2-D model of the sky for optimal extraction
  *	ap_trace():  Traces spectrum automatically or manually / interactively
  *	optimal():  Optimal extraction
  *	ap_extract():  Summed extraction, sky subtraction, generates S/N plot
  *	wavelength():  Apply wavelength calibration
  *	flux_calibration():  Applies flux calibration

#### <ins>spec.py:</ins>

This script defines all the input parameters and imports the functions from Spectroscopic_routines.py.  This script can be setup for a particular setup / instrument.  See example script.

Currently all the scripts and files have to be in the same folder as your spectra.  **This will be changed eventually**

This is the script you run:

***python spec.py***

### 4. Required data files:
  * <ins>Reference spectrum</ins> (wavelength and counts columns).  The script assumes angstrom units.  See ref_cuar_gr7_162.dat as example.
  * <ins>Line list</ins> (wavelength and line ID).  The script assumes angstrom units.  See line_list.dat as example.

