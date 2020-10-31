%README.txt

%{
1. Begin by running 

>> generateDeconvolutionFilesForMidgetRGCs(args)

which computes how:
- the RF center (constructed of 1-N cones, each with a weight of 1),and
- the RF surround (constructed of M cones with Gaussian weights)
radius and peak sensitivity are changed by applying a subject's PSF at different 
eccentricities. The change in sigma and amplitude provide the data for the deconvolution
model which allow us to transform the Croner&Kaplan '94 visual RF parameters 
(radius and sensitivity) from visual to retinal counterparts and vice-versa.
See generateDeconvolutionFilesForMidgetRGCs() for different ways to run this function

2. Once a set of deconvolution files for the center and the surround subregion
have been generated, say for eccentricities [0 1 2 3 4],  we need to 
- either update the constant property @CronerKaplanRGCModel.defaultDeconvolutionEccs 
with that eccentricity range, or 
- instantiate @CronerKaplanRGCModel objects by passing the ('deconvolutionEccs', [0 1 2 3 4]) key-value pair.

%}
