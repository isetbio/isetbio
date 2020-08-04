%README.txt

%{ 
We begin by running 

>> performPSFConvoComputations('generateNewDeconvolutionFiles', true)

which convolves a number of 2D Gaussian shapes with the PSFs at different eccentricities
and quantifies the change in spread and amplitude reduction of the input Gaussian. 
This information simulates how a retinal Gaussian subregion (center/surround)
would be changed in visual space (say when mapped in an in-vivo preparation).
The change in sigma and amplitude provide the data for the deconvolution
model which allow us to transform the Croner&Kaplan '94 visual RF parameters 
(radius and sensitivity) from visual to retinal counterparts and vice-versa.

The generated data are saved in the CronerKaplanRGCModel.psfDeconvolutionDir,
with a separate file for each of the eccentricities examined.

To visualize an existing dataset, call:

%}
