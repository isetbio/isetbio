%README.txt

%{ 
We begin by running 

>> performPSFConvoComputations('generateNewDeconvolutionFiles', true)

which computes how:
- the RF center (constructed of 1-N cones, each with a weight of 1),and
- the RF surround (constructed of M cones with Gaussian weights)
radius and peak sensitivity are changed by applying a subject's PSF at different 
eccentricities. The change in sigma and amplitude provide the data for the deconvolution
model which allow us to transform the Croner&Kaplan '94 visual RF parameters 
(radius and sensitivity) from visual to retinal counterparts and vice-versa.

The generated datasets are saved in the CronerKaplanRGCModel.psfDeconvolutionDir,
with a separate file for each of the eccentricities examined.

To visualize an existing dataset, for a subset of the eccentricities, call:
performPSFConvoComputations('eccTested', [0 0.2 0.5 1 1.5 2 3 4])

%}
