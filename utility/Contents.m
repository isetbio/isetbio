% UTILITY
%
% Top Level Files
%   checkToolbox         - Checks whether certain matlab toolbox has been installed
%   dpi2mperdot          - Convert dots per inch to meters format
%   getMiddleMatrix      - Extract values near middle of a matrix.
%   hist2d               - Calculate and return a 2D histogram of D
%   ieClip               - Clip data to range specified arguments.
%   ieCmap               - Prepare simple color maps
%   ieContrast           - Convert a signal to a contrast represenation
%   ieDeg2rad            - Convert degrees to radians
%   ieFindWaveIndex      - Returns a (0/1) vector of indices s.t. wave(idx) matches a waveVal entry.
%   ieGIF                - Save movie data (x, y, t) as a gray scale GIF
%   ieHwhm2SD            - Convert half width half max to standard deviation for Gaussian
%   ieHwrect             - Half-wave rectification of the signal data
%   ieMovie              - Show a movie of an (x, y, t) or (x, y, c, t) matrix
%   ieN2MegaPixel        - Compute megapixel count from N
%   ieRad2deg            - Convert radians to degrees
%   ieScale              - Scale the value in im into a specified range
%   ieShape              - Draw a shape on the current window
%   ieSpace2Amp          - Transform spatial data to amplitudes in cycles per spatial unit  
%   ieSpaceTimeFilter    - Apply a filter kernel to the space-time signal
%   ieUncompressData     - Deprecated.
%   ieUnitScaleFactor    - Return scale factor that converts from meters or seconds to other scales
%   ieWave2Index         - Convert a wavelength to an index into the wave list.
%   isetbioDataPath      - Return path to oldstyle directory containing data bundled with isetbio
%   isetbioRootPath      - Return the path to the root isetbio directory
%   isetRootPath         - Return the path to the isetbio isettools directory
%   isodd                - Always nice to know if there is something odd going on ;).
%   notDefined           - Test whether a variable (usually a function argument) is defined
%   sample2space         - Return the physical position of samples
%   temporalEquivEcc     - (IN PROG.) Equivalent temporal eccentricity acc. to Kalmar/Chichilnisky
%   twoGammaResp         - Create temporal response composed of the difference of 2 gamma functions
%   unitFrequencyList    - Calculate a vector of normalized frequencies for an N-vector
%   unitLength           - Convert vectors in the rows of v to unit length in rows of u
%   unpadarray           - Inversion of padarray
%   upperQuad2FullMatrix - Duplicates the upper right quad of a matrix in the other quads
%
% FILE (subdir file)                    File utility functions for isetbio
%   file/ieImageType                  - Determine the type of image in a file
%   file/ieReadSpectra                - Read in spectral data and interpolate to the specified wavelengths
%   file/ieSaveMultiSpectralImage     - Save a Matlab data file containing a multi-spectral image.
%   file/ieSaveSIDataFile             - Write file with data for shift-invariant optics
%   file/ieSaveSpectralFile           - Save a spectral data ISET data file
%   file/ieVarInFile                  - Check whether a variable is in a particular Matlab file
%   file/vcExportObject               - Save an object into the relevant object directory
%   file/vcImportObject               - Import an ISET structure from a file to the vcSESSION structure 
%   file/vcReadImage                  - Read image color data, return multispectral photons
%   file/vcSaveObject                 - Deprecated. Use vcExportObject. Save an ISET object into a .mat file. 
%   file/vcSelectDataFile             - Select a data file name for reading or writing
%   file/vcSelectImage                - Return the full path name and type of an ISET image.
%
% FUNCTIONS (subdir functions)          Computational utility functions for isetbio
%   functions/coneEmpiricalDimFlash   - Return estimate of Schnapf et al. dim flash cone photocurrent response.
%   functions/hill                    - Compute hill function
%   functions/sumOfTwoDblExponentials - Compute sum of two double exponential functions, with no mean offset.
%   functions/weberFechner            - Calculate the Weber Fechner tvi function
%
% HYPERCUBE (subdir hypercube)          Hypercube-related functions in isetbio
%   hypercube/hcBasis                 - Create wavelength basis functions and coefficients for an hc image
%   hypercube/hcBlur                  - Blur each plane in a hypercube with conv2
%   hypercube/hcIlluminantScale       - Estimate the relative illuminant intensity across space
%   hypercube/hcimage                 - Display a hypercube image.
%   hypercube/hcimageCrop             - Select region to crop from a hypercube image
%   hypercube/hcimageRotateClip       - Clip and rotate hypercube data. Used for visualization with specularities
%   hypercube/hcReadHyspex            - Reads an ENVI image.
%   hypercube/hcReadHyspexImginfo     - Read ENVI image header files (special purpose for hyperspectral data)
%
% IMAGE (subdir image)                  Image operations for isetbio
%   image/convolvecirc                - Performs 2D circular convolution 
%   image/ieLUTDigital                - Convert DAC values to linear RGB values through a gamma table
%   image/ieLUTInvert                 - Calculate inverse lookup table (lut) at certain sampling steps
%   image/ieLUTLinear                 - Convert linear RGB values through a gamma table to DAC values
%   image/imageContrast               - Compute the image contrast in an RGB style image
%   image/imageFlip                   - Flip image data   - updown or leftright 
%   image/imageHarmonic               - Creates a sum of harmonic images, potentially windowed by a Gaussian
%   image/imageIncreaseImageRGBSize   - Increase the size of an rgb-style image (r, c, w) by pixel replication.
%   image/imageLinearTransform        - Apply a linear transformation to the channels of an RGB or XW image 
%   image/imageMakeMontage            - Create image data comprising a montage of the slices in hc
%   image/imageMontage                - Create a window with a montage of the slices in the hypercube data
%   image/imageMultiview              - Display multiple images of selected GUI objects
%   image/imageRotate                 - Rotate image data   - CW or CCW
%   image/imagescRGB                  - Display a scaled RGB format image. 
%   image/imageSPD                    - Derive an RGB image from an SPD image
%   image/imageTranspose              - Transpose image data in multispectral or RGB
%
% PLOTS (subdir plots)                  Various isetbio plots
%   plots/identityLine                - Draw an identity line on the current axis
%   plots/ieFormatFigure              - Formats a figure for presentations or papers.
%   plots/oiPlot                      - Gateway routine for plotting optical image (oi) properties
%   plots/plotSetUpWindow             - Initialize position, size and colormap for ISET plot window
%   plots/plotSpectrumLocus           - Draw the outline of spectrum locus on the chromacity diagram
%   plots/plotTextString              - Add a text string to a 2D graph in a specific position.
%   plots/scenePlot                   - Gateway routine to plot scene radiance properties
%
% PROGRAMMING (subdir programming)      Programming operations for isetbio
%   programming/addText               - Add text to an existing text string
%   programming/cellDelete            - Delete some cells from a list of cells 
%   programming/checkfields           - Check for the presence of a field, possibly nested, in a structure
%   programming/clx                   - Script to clear many things
%   programming/gatherStruct          - Gather distributed struct to current working directory
%   programming/hiddenHandle          - Empty abstract class inheritated from handle class
%   programming/ieInputParser         - Get an inputParser configured with ISETBIO conventions
%   programming/ieMemoryCheck         - Calculate the memory use in the calling routine or base
%   programming/ieParameterOtype      - Determine object type from parameter name
%   programming/ieParamFormat         - Converts s to a standard ISET parameter format  (lower case, no spaces)
%   programming/ieRemovePathDir       - Removes any .svn directories from the pathList.
%
% PSYCHOPHYSICS (subdir psychophysics)  Related to calculating psychophysical performance. The TSD
%                                       routines are based on the equal-variance normal assumption.
%   psychophysics/analyticpHitpFa     - Compute hit and fa rates from underlying distributions and criterion
%   psychophysics/analyticPoissonIdealObserver - Performance of analytic ideal observer
%   psychophysics/computeDPrimeCritNorm - Compute d-Prime and criterion from hit and fa rates
%   psychophysics/computeROCArea      - Compute area under ROC curve from underlying distributions
%   psychophysics/dPrimeToTAFCFractionCorrect - Calculate TAFC fraction correct from d-prime
%   psychophysics/TAFCFractionCorrectToDPrime - Calculate d-prime from TAFC fraction correct
%
% STATISTICS (subdir statistics)        Statistical calculations for isetbio
%   statistics/biNormal               - Compute bivariate normal function
%   statistics/expRand                - Exponentially distributed random number generator
%   statistics/gammaPDF               - Gamma function   - used for impulse response calculations and HIRF
%   statistics/getGaussian            - Create a 2D Gaussian distribution pdf
%   statistics/ieExprnd               - Random arrays from exponential distribution.
%   statistics/ieMvnrnd               - Multivariate normal random number generator
%   statistics/ieNormpdf              - Normal probability density function (pdf).
%   statistics/iePoisson              - Create a matrix of Poisson samples using rate parameters in lambda
%   statistics/iePrctile              - Measures the percentiles of the sample in X.
%   statistics/lorentzSum             - Compute sum of Lorentzian components
%


