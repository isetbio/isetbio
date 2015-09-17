Notes on mosaic code


The Curcio/Sloan model is for the overall cone packing, without regard to cone type.  It is implemented in the function MakeSloanList, which takes an argument for the number of pixels (given as linear number of pixels for a square image), the number of cones, and the parameters of the model.

The model and how parameters vary with eccentricity is described in: Curcio CA, Sloan KR. 1992. Packing geometry of human cone photoreceptors: variation with eccentricity and evidence for local anisotropy. Visual Neurosci 9: 169-80.

To generate a trichromatic mosaic, the routine is MakeSloanSClist.  This takes a total number of cones, plus the number that should be S and M.  It generates an overall mosaic, then a mosaic with the same parameters with the number of S cones.  The cones in the overall list that are nearest to the generated S cones are called S.  The rest are allocated between L and M by coin flipping according to the L:M ratio implied by the input parameters.

The output of MakeSloanList is a list of x,y positions.  The output of MakeSloanClist is x,y positions plus an indicator variable (1 -> L, 2 -> M, 3 -> S) giving type of each cone.

This code was written in the early 1990's and was part of a larger toolbox that did various image manipulations related to sampling.  But I think ISETBIO handles most of those, so I have just tried to strip out the key stuff for creating the mosaics.  Matlab had only two-dimensional matrices at the time this was written, so some things may be coded in what seems like a convoluted way given modern Matlab.

