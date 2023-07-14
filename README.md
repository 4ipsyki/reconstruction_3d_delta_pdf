Procedures and algorithms for 3D-(delta)PDF data analysys are briefly described.
Scripts as ip.m, find_ub.m, calc_dk_P.m and imga2hkl.m are written by
M. v. Zimmermann (DESY), and esthetically modified by O. Ivashko (DESY).
The rest of the scripts are written by O. Ivashko (DESY).

July 2023
================================================================================
The following list is orginized in an ordered fashion and serves only as a
generalized guideline. Scripts and relative parameters need to be adjusted
case-by-case, as most of the algorithms are not automatically self-tuned.
These scripts are polished from analyses on several data taken with PILATUSX 
CdTe 2M and Perkin-Elmer XRD1621 at P07-EH2 and P21.1 at PETRAIII, DESY.
The general aquisition mode consists in several detector macro-positions,
in order to cover the largest Q-range possible, and three or four micro-positions,
to cover the gaps of the detector (when present). The sample is rotated continuosly
by 180+ or 360+ degrees while firing the detector with a chosen exposure time.
The motor speed, and thus the delta-angle, is adjusted as needed. Experiments are
performed in a forward transmission geometry giving the high energy of X-rays of
the order of 100 keV. To reduce the background, a scattering chamber with internal
beamstop and filled with He or completely evacuated is typically used. Giving an
opening angle of ~ +- 30 gedrees and 100 keV, total scattering measurements are
possible since the Q range extends to ~ 27 \AA^-1.
================================================================================

<> calibrate detector parameters.
	Usually a powder capillary filled with LaB6, CeO2 or Ni is measured at the
	sample position. When using pyFAI-calib2, a convertertion to beam-origin
	coordinates must be performed.
	--> use poni2org.m to convert from the poni file.

<> create [input_parameters].mat file
	This contains the neccesarry parameters and directories for the reconstruction.
	--> use ip.m as an example.

<> (OPTIONAL) create custom mask.
	This must contain bad pixels and gaps of the detector.
	Sometimes also a considerable portion of the border must be masked.
	Include beamstop and its holder as well as the border of the chamber.
	--> use create_mask.m as an example.
		NOTE: a basic mask (dead pixels / blank) must be incuded anyway.

<> generate UB (orientation) matrix.
	This is a matrix containing the coordinates of the crystal axes (in real space)
	stacked in collumns and referred to nominal angle posiitons. The convention 
	is the same as in XDS, and it's matrix can be directly used.
	--> use find_ub.m as an example.
		--> matrix_fun.m
		--> matrix_constrain.m
		(--> matrix_fun_detector.m      )
		(--> matrix_constrain_detector.m)

<> (OPTIONAL) extract/create absorption correction.
	This includes in a general way all the isotropic frame-by-frame corrections.
	It can include the variation of the incident beam or the change of volume
	fraction of the sample in the beam.
	--> use absorption_correction.m as an example.
		NOTE: when a direct-beam monitor after the sample is present,
			  the above-mentioned script is obsolete.

<> perform the reconstruction.
	Here the previously choosed and determined parameters and corrections have to
	be included. Partial or complete reconstructions can be performed and saved.
	--> use reconstruct.m as an example.
		--> calc_dk_P.m
		--> image2hkl.m
		(--> merge.m)

<> (OPTIONAL) plot the reconstruction.
	Basic plotting along the three principle axes with choosen Q point and
	out-of-plane binning.
	--> use fast_plot.m as an example.
		--> slice.m
		--> savepng.m

<> (OPTIONAL) perform symmetrization.
    Application of Laue symmetries on the dataset. This step helps to reduce the amount
    of missing data, as due to assymetric detector position, beamstop holder, etc.
    This step can be also done after the punch & fill part (for 3D-deltaPDF).
    --> use perform_symmetry.m as an example.
        --> apply_symmetry.m
        --> compare_plot.m

########### ======= FOR 3D-deltaPDF ======= ###########

<> perform Bragg-peak punching.
	A statistical method for outlayer detection is adopted (Ref. KAREN).
	Within the algorith, also a statistical and fast filling can be performed.
	--> use perform_punch.m as an example
		--> paf_karen.m
        --> compare_plot.m
        
<> perform filling.
	A Gaussian-core convolution filter is used based on the python Astropy 	
	package (Ref. Astropy). This method is chosen for simplicity, speed and
	overall consistency of the rersult. A matlab script is used for comodity,
	but a python version is easily implementable.
	--> use perform_fill.m
		--> fill_astropy.m
        --> compare_plot.m
    
(<> (OPTIONAL) perform symmetrization.)
    
<> perform iFFT
	
	--> use perform_ifft.m

<> (OPTIONAL) plot the 3D-deltaPDF

	--> use fast_plot_pdf.m as an example.
		--> slice.m
		--> savepng.m






