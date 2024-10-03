The Image System Engineering Toolbox for Biology - ISETBio - is a Matlab toolbox for calculating the properties of the front end of the visual system.  The toolbox begins with a description of the scene radiance, how the radiance is transformed by the optics into the retinal image, captured of light by the photopigment, the photocurrent responses in the receptors.  We are actively working on modeling how the photoreceptor signals are converted into retinal ganglion cell responses.

This repository includes a [WIKI](https://github.com/isetbio/isetbio/wiki) that describes the software as well as many examples of how to perform computations to calculate the visual encoding of light in the eye.  The [WIKI](https://github.com/isetbio/isetbio/wiki) also describes tools to quantify, view and and analyze the information contained at different neural stages.  

As of May 27, 2024 To run ISETBio, you must download ISETCam and have it on your Matlab path.  See installing [ISETCam](https://github.com/iset/isetcam/wiki). 

A set of tutorial videos introducing aspects of ISETBio and ISETCam is available at [this YouTube playlist](https://www.youtube.com/playlist?list=PLr6PuubdQrtQ-rz5RIe13k3YFrmwBh7gr).

### History 

* May 29, 2024 - See changes to ISETCam in the ReadMe for that repository. The validation routines for ISETBio have been moved into a new [ISETValidations](https://github.com/ISET/isetvalidate) repository.

### Ancient history

The ISETBIO code includes a portion of Image Systems Engineering Toolbox (ISET) that is sold by [Imageval Consulting, LLC](http://www.imageval.com).  That code is designed to help industrial partners design novel image sensors. The ISETBIO portion of the ISET code is freely distributed for use in modeling image formation in biological systems. 

ISETBIO also includes the WavefrontOptics code developed by David Brainard, Heidi Hofer and Brian Wandell.  That code implements methods for taking adaptive optics data from wavefront sensors and calculating the optical blur as a function of wavelength for model human eyes.  The toolbox relies on data collected by Thibos and colleagues. We also gratefully acknowledge important contributions from Jon Winawer.  
