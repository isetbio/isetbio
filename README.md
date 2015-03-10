The Image System Engineering Toolbox for Biology ISETBIO is a Matlab toolbox designed for calculating the properties of the front end of the visual system.  The ISETBIO code includes a portion of Image Systems Engineering Toolbox (ISET) that is sold by Imageval Consulting, LLC.  That code is designed to help industrial partners design novel image sensors. The ISETBIO portion of the ISET code is freely distributed for use in modeling image formation in biological systems. 

ISETBIO also includes the WavefrontOptics code developed by David Brainard, Heidi Hofer, and Brian Wandell.  That code implements methods for taking adaptive optics data from wavefront sensors and calculating the optical blur as a function of wavelength for model human eyes.  The toolbox relies on data collected by Thibos and colleagues.

ISETBIO also includes a set of tools, cones and rgc, implemented by Jon Winawer and Brian Wandell to transform cone absorptions into retinal ganglion cell spikes.  

The ISET and WavefrontOptics portions are reasonably mature. The cones and rgc code are in the very early development stage.

For a general introduction to human vision, please see:

  http://foundationsofvision.stanford.edu/

For examples of ISET/VSET code and tutorials on human vision, please see:

  https://www.stanford.edu/group/vista/cgi-bin/FOV/computational-examples/

        Copyright (c) 2012, Brian A. Wandell wandell@stanford.edu

Implementation notes

The toolbox is written around several data structures.  Each is implemented using set/get/create syntax.  Major computations are implemented by computing with these data structures.

The scene data structure describes the scene radiance (photons).  It is set up to permit depth encoding, though nearly all of the current examples are based on a radiance field originating from a single plane.

The optical image transforms the irradiance distribution at the sensor, after the radiance has passed through the topics.  There are several computational models that implement the transformation: diffraction limited, shift invariant, and ray trace (shift variant). You can choose the one you want to use depending on the level of information you have about the optics.  There is an implementation of the optics of the human eye, based on work from Marimont and Wandell.  

This toolbox is coordinated with the Wavefront Toolbox (also in github) uses data from adaptive optics to simulate the defocus based on measurements of many different human eyes.

The sensor transforms the irradiance into a spatial array of cone absorptions. The sensor pixels and spectral quantum efficiency can be set to model the human cones at various eccentricities and with various types of inert pigments (macular pigment, lens density, optical density).  The sensor simulates the photon absorptions in the cone (and rod) mosaics.

