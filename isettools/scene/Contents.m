% SCENE
%
% Top Level Files
%   sceneAdd                - Add together the photons from two scenes
%   sceneAdjustIlluminant   - Adjust the current scene illuminant to the value in data
%   sceneAdjustLuminance    - Scale scene mean luminance to meanL
%   sceneAdjustReflectance  - Scale the data so that the peak reflectance is 1
%   sceneCalculateLuminance - Calculate scene luminance (cd/m^2)
%   sceneClearData          - Clear the scene data entries. 
%   sceneComplete           - Complete the scene parameters (geometry and spatial) after initialization
%   sceneCreate             - Create a scene structure.
%   sceneCrop               - Crop scene data.
%   sceneDescription        - Text description of the scene properties, displayed in scene window
%   sceneExtractWaveband    - Extract wave bands from the scene
%   sceneFrequencySupport   - Compute spatial frequency support for scene in specified units 
%   sceneFromBasis          - Create a scene from a structure with basis functions (linear model)
%   sceneFromFile           - Create a scene structure by reading data from a file
%   sceneGet                - Get scene parameters and derived properties
%   sceneInit               - First stage in building up a scene
%   sceneInitGeometry       - Initialize scene distance parameter.  
%   sceneInitSpatial        - Initialize the scene field of view to 10 deg.  
%   sceneInterpolateW       - Wavelength interpolation for scene image data
%   sceneList               - Scene types with parameters
%   sceneSaveImage          - Write png image approximating appearance of photon data
%   sceneSet                - Set ISET scene parameter values
%   sceneShowImage          - Render an image of the scene data
%   sceneSpatialResample    - Spatial resample all wavebands of a scene
%   sceneSpatialSupport     - Calculate the spatial positions of the scene sample locations
%   sceneSPDScale           - Multiply, Divide, Add or Subtract the scene radiance data
%   sceneToFile             - Write scene data in the hyperspectral and multispectralfile format
%   sceneWBCreate           - Create a directory of waveband scene images in separate files  
%
% ADATPIVEOPTICS (subdir adaptiveOptics)
%   AOMonochromaticCornealPowerToRadiance - Convert power at cornea to
%   equivalent radiance.
%
% ILLUMINATION (subdir illumination) isetbio scene illumination functions
%   illumination/blackbody           - Generate the spectral power distribution of a blackbody radiator
%   illumination/illuminantCreate    - Create an illuminant (light source) structure. 
%   illumination/illuminantGet       - Get parameter value from an illuminant structure
%   illumination/illuminantModernize - Convert old format illuminant structures to the modern format
%   illumination/illuminantRead      - Return spectral radiance of a standard illuminants in energy units
%   illumination/illuminantSet       - Set parameter value for illuminant structure
%
% IMGTARGETS (subdir imgtargets)       Image Targets
%   imgtargets/imgMackay             - Create a MacKay chart spatial pattern. 
%   imgtargets/imgRadialRamp         - Make a radial ramp function
%   imgtargets/imgRamp               - Create a set of intensity ramps as test spatial pattern. 
%   imgtargets/imgSweep              - Create a sweep frequency image as a test pattern.
%   imgtargets/imgZonePlate          - Make a zone plate image
%
% MACBETH (subdir macbeth)             Macbeth Color Chart scene structures
%   macbeth/macbethChartCreate       - Initiate a scene structure of the Gretag/Macbeth Color Chart reflectances
%   macbeth/macbethColorError        - Analyze color error of a MCC in the image processor window
%   macbeth/macbethCompareIdeal      - Create an image of an ideal MCC (color, temp, ...) with data embedded
%   macbeth/macbethEvaluationGraphs  - Evaluate linear fit L from sensor rgb to the ideal rgb of an MCC
%   macbeth/macbethGretagSGCreate    - (Non-function) Create a scene data structure of the Gretag digital color chart SG
%   macbeth/macbethIdealColor        - Calculate MCC values for given an illuminant in a color space
%   macbeth/macbethLuminanceNoise    - Analyze luminance noise in gray series of MCC from image processor window
%   macbeth/macbethPatchData         - Return a cell array with the linear RGB values from a vcimage or sensor 
%   macbeth/macbethReadReflectance   - Read the macbeth surface reflectances into the standard vimage ordering
%   macbeth/macbethRectangles        - Calculate mid point of rectangles for an MCC from the corner points
%   macbeth/macbethROIs              - Derive a rectangular region of interest for an MCC chart
%   macbeth/macbethSensorValues      - Identify MCC patches and calculate RGB mean (s.d.) of the 24 patches
%
% PATTERN (subdir pattern)             Scene Patterns
%   pattern/FOTarget                 - Create a frequency/orientation image
%   pattern/harmonicP                - Default harmonic params
%   pattern/ieCheckerboard           - Create a checkerboard image
%   pattern/sceneGridLines           - Create scene comprising an array of grid lines
%   pattern/sceneHarmonic            - Create a scene of a Gaussian windowed harmonic function (Gabor patch)
%   pattern/scenePointArray          - Make a point array stimulus for evaluating the optics
%   pattern/sceneRamp                - Intensity ramp (see L-star chart for L* steps)
%   pattern/sceneReflectanceChart    - Create a reflectance chart for testing
%   pattern/sceneSlantedBar          - Set the scene to equal photons across the wavelength(s)
%   pattern/sceneVernier             - scene for vernier acuity
%   pattern/vernierP                 - Default vernier stimulus parameters
%
% REFLECTANCE (subdir reflectance)     Scene reflectances
%   reflectance/ieReflectanceSamples - Return a sample of reflectances 
%
% SCENEDEPTH (subdir scenedepth)       Scene depth functions
%   scenedepth/sceneDepthRange       - Create a scene with photons restricted to a particular depth range
%
% SCENEGUI (subdir scenegui)           Scene GUI functions
%   scenegui/sceneClose              - sceneClose - close SCENE window.
%   scenegui/sceneInterpolate        - Spatially interpolate the scene radiance data by sFactor
%   scenegui/sceneKeyPress           - Determine scene window key press callbacks.
%   scenegui/sceneSetEditsAndButtons - Fill scene window fields based on the current scene information
%   scenegui/sceneSetRowCol          - GUI to set scene row and col (image size)
%   scenegui/sceneWindow             - Graphical user interface to manage the ISET SCENE properties.
%
