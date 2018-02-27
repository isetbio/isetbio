% CONES
%
% Top level
%   coneMosaicWindow - Cone image coneMosaicWindow interface
%
% Cone Mosaic (subdir @coneMosaic)
%   @coneMosaic/applyEMPath                        - Apply eye movement path and return absorptoins absorptions.
%   @coneMosaic/clearData                          - Clear computed data
%   @coneMosaic/compute                            - compute  Method of coneMosaic, to compute the cone absorptions, possibly for multiple trials (repeats)
%   @coneMosaic/computeCurrent                     - Convert absorptions to photocurrent using the os model.
%   @coneMosaic/computeForOISequence               - Compute cone absorptions and optionally photocurrents for a @oiSequence
%   @coneMosaic/computeSingleFrame                 - Single frame compute function for coneMosaic object.
%   @coneMosaic/coneAbsorptions                    - Get isomerizations from one of the cone classes
%   @coneMosaic/coneMosaic                         - coneMosaic  Create a cone mosaic object
%   @coneMosaic/conePositions                      - Return the cone sample positions in meters
%   @coneMosaic/copyElement                        - Make a shallow copy of the properties
%   @coneMosaic/demosaicedResponses                - Return demosaiced response maps (absorptions or currents) from a coneMosaic object
%   @coneMosaic/description                        - Create string with key mosaic parameters
%   @coneMosaic/descriptionOS                      - DESCRIPTIONSOS  Summarize properties of the outer segment
%   @coneMosaic/emGenSequence                      - Generate sequence of eye movements
%   @coneMosaic/lowPassMosaicResponse              - lowPassMosaicResponse  Low pass filter mosaic isomerizations
%   @coneMosaic/photonNoise                        - Add photon noise to the absorptions 
%   @coneMosaic/plot                               - Plot function for @conemmsaic base class
%   @coneMosaic/retinalCoverage                    - Compute aperture and geometric retinal coverage
%   @coneMosaic/setSizeToFOV                       - Updates the cone mosaic size to match the FOV
%   @coneMosaic/setWave                            - Keep photopigment and macular pigment wavelength sampling consistent.
%   @coneMosaic/tResample                          - Resample an absorption time sequence
%   @coneMosaic/window                             - Opens a cone mosaic window GUI.
%
% Cone Mosaic Hex (subdir @coneMosaicHex)
%   @coneMosaicHex/computeActivationDensityMap     - Compute activation images for the hex mosaic (all cones + LMS submosaics)
%   @coneMosaicHex/computeDensityMap               - Method to compute the cone density of @coneMosaicHex
%   @coneMosaicHex/coneMosaicHex                   - Create a hexagonal cone mosaic class
%   @coneMosaicHex/description                     - Object description for the cone mosaic hex case
%   @coneMosaicHex/displayInfo                     - Print various infos about the cone mosaic
%   @coneMosaicHex/movieHex                        - Generates a movie of the coneMosaicHex absorptions or current over time.
%   @coneMosaicHex/plot                            - Plot function for coneMosaicHex (subclass of coneMosaic)
%   @coneMosaicHex/plotHexMeanImage                - Plot mean absorptions or current from cone mosaic hex into a window image
%   @coneMosaicHex/plotHexMosaic                   - Visualize different aspects of the hex grid
%   @coneMosaicHex/plotMosaicProgression           - Plot the mosaic progression
%   @coneMosaicHex/reassignConeIdentities          - Reassign the cone identities of the cone mosaic hex object
%   @coneMosaicHex/renderActivationMap             - Render (in the passed axesHandle) an activation map for the hex mosaic
%   @coneMosaicHex/renderHexMesh                   - Render the hex mesh for the cone mosaic hex object
%   @coneMosaicHex/renderPatchArray                - Render the patch array
%   @coneMosaicHex/resampleGrid                    - Sample the original rectangular mosaic using a hex grid sampled at the passed resamplingFactor
%   @coneMosaicHex/reshapeHex2DmapToHex3Dmap       - Reshape the existing 2D Hex Map to a 3D Hex Map
%   @coneMosaicHex/reshapeHex3DmapToHex2Dmap       - Reshape the existing 3D Hex Map to a 2D Hex Map
%   @coneMosaicHex/restoreOriginalResState         - Restore the original Rectangular Grid State
%   @coneMosaicHex/saveOriginalResState            - Save the original Rectangular Grid State
%   @coneMosaicHex/setSizeToFOVForHexMosaic        - Set the size to the field of view for the hex mosaic
%   @coneMosaicHex/visualizeActivationMaps         - Separately visualize mosaic activations for each submosaic and the whole
%   @coneMosaicHex/visualizeGrid                   - Visualize different aspects of the hex grid
%
% Macular (subdir @Macular)
%   @Macular/description                           - Generate description string for this macular pigment
%   @Macular/Macular                               - Class for macular pigment properties
%
% Photo Pigment (subdir @photoPigment)
%   @photoPigment/description                      - Photopigment text description
%   @photoPigment/photoPigment                     - Class for single cone photopigment and related properties
%
% Outer Segment (subdir outersegment)
%   Utilities (subdir utilities)
%       outersegment/utilities/osCreate            - Create an outer segment object
%       outersegment/utilities/osAddNoise          - Add noise to noise-free membrane current
%
%   Outer Segment (subdir @outerSegment)
%       outersegment/@outerSegment/osGet           - Get base class outer segment parameters
%       outersegment/@outerSegment/osSet           - Sets isetbio outersegment object properties for base class
%       outersegment/@outerSegment/outerSegment    - Parent class for computing cone photocurrent (pA) from isomerization rate (R*)
%       outersegment/@outerSegment/plot            - Plot for the outersegment base class
%       outersegment/@outerSegment/resample        - Method to resample the photocurrents using bicubic interpolation.
%
%   Outer Segment Linear (subdir @osLinear)
%       outersegment/@osLinear/linearFilters       - Returns the photocurrent impulse response for a single absorption
%       outersegment/@osLinear/osCompute           - Linear model computing outer seg. photocurrent from isomerizations (R*) 
%       outersegment/@osLinear/osGet               - Gets isetbio outersegment object parameters.
%       outersegment/@osLinear/osLinear            - Linear subclass of OS object, convert isomerizations (R*) to current
%       outersegment/@osLinear/osSet               - Sets isetbio outersegment object parameters.
%       outersegment/@osLinear/plot                - Plot osLinear properties or sends to @outersegment.plot
%
%   Outer Segment Identity (subdir @osIdentity )
%       outersegment/@osIdentity/osGet             - Gets the isetbio outersegment object parameters.
%       outersegment/@osIdentity/osIdentity        - A subclass of @outerSegment object
%       outersegment/@osIdentity/osSet             - Sets isetbio outersegment object parameters.
%
%   Outer Segment Display RGB (subdir @osDisplayRGB)
%       outersegment/@osDisplayRGB/osCompute       - Pass the cone isomerizations (R*) without any temporal filtering. 
%       outersegment/@osDisplayRGB/osDisplayRGB    - displayRGB subclass of the outersegment object
%       outersegment/@osDisplayRGB/osGet           - Gets the isetbio outersegment object parameters.
%       outersegment/@osDisplayRGB/osSet           - @osDisplayRGB that sets OS object parameters using the input parser.
%
%   Outer Segment BioPhysical (subdir @osBioPhys)
%       outersegment/@osBioPhys/osAdaptSteadyState - Steady-state background current calculated from the background rates.
%       outersegment/@osBioPhys/osAdaptTemporal    - Time varying current response from photon rate and initial state
%       outersegment/@osBioPhys/osBioPhys          - Create a biophysically based outersegment (os) object.
%       outersegment/@osBioPhys/osCompute          - Compute the photocurrent response of the L, M and S cones
%       outersegment/@osBioPhys/osGet              - Gets the isetbio outersegment object parameters.
%       outersegment/@osBioPhys/osSet              - Sets the isetbio outersegment object parameters.
%
% Utilities (subdir utilities)
%   utilities/coneMeanIsomerizations                                - Calculate spatial mean photon rate(R*/sec) for the 3 cone types in mosaic
%   utilities/coneTypeLocations                                     - Return row/col or indices of the three types of cones
%   utilities/diameterForCircularApertureFromWidthForSquareAperture - Calculate diameter of circular aperture of identical area given width for
%   utilities/sizeForSquareApertureFromDiameterForCircularAperture  - calculate width of square aperture of identical area given diameter for
%
