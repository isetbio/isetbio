function t_mRGCMosaicVisualize(options)
% Introduction to the midget RGC mosaic (mRGCMosaic) object.
%

% Requires locally configured precomputed mRGC resources.
% SkipFile
% Description:
%    Demonstrates
%        - how to visualize different aspects
%


% History:
%    10/25/23  NPC  Wrote it.
%    7/26/2026 NPC  Modernized it

% Examples:
%{

    t_mRGCMosaicVisualize();

    t_mRGCMosaicVisualize(...
        'rgcMosaicName', 'JCNpaperTemporal7DegsMosaic', ...
        'opticsSubjectName', 'JCNpaperDefaultSubject');

%}


arguments

    % ---- Mosaic specifiers for selecting a prebaked mRGC mosaic ---

    % See RGCMosaicConstructor.helper.utils.initializeRGCMosaicGenerationParameters
    % for what is available and to add new mosaics
    options.rgcMosaicName (1,:) char = 'JCNpaperTemporal4DegsMosaic';


    % ---- Which species to employ ----
    % Choose between {'macaque', 'human'}. If 'macaque' is chosen, the input
    % cone mosaic has a 1:1 L/M cone ratio.
    options.coneMosaicSpecies  (1,:) char {mustBeMember(options.coneMosaicSpecies,{'human','macaque'})} = 'human';

     % ----- Which subject optics to employ -----
    options.opticsSubjectName (1,:) ...
        char ...
        {...
        mustBeMember(options.opticsSubjectName, ...
            { ...
            'JCNpaperDefaultSubject' ...
            'JCNpaperSecondSubject' ...
            'VSS2024TalkFirstSubject' ...
            'VSS2024TalkSecondSubject' ...
            'JCNpaperStrehlRatio_0.87' ...
            'JCNpaperStrehlRatio_0.72' ...
            'JCNpaperStrehlRatio_0.59' ...
            'JCNpaperStrehlRatio_0.60' ...
            'JCNpaperStrehlRatio_0.27' ...
            'JCNpaperStrehlRatio_0.23' ...
            'JCNpaperStrehlRatio_0.21' ...
            'JCNpaperStrehlRatio_0.19' ...
            'JCNpaperStrehlRatio_0.09' ...
            } ...
            ) ...
        } ...
        = 'JCNpaperDefaultSubject';


    % ------ targetVisualSTF options ----
    % Options are : {'default', 'x1.3 RsRcRatio'}
    % These are with respect to the macaque data of the Croner & Kaplan '95 study
    % 'default': target the mean Rs/Rc, and the mean Ks/Kc (Rs/Rc)^2
    % See RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.generateTargetVisualSTFmodifiersStruct
    % for all existing options
    options.targetVisualSTFdescriptor (1,:) char = 'default';


    % ------ Visualization options ----
    % Visualize cone pooling maps for a target RGC
    options.targetRGCindexForVisualizingConePoolingMaps (1,:) double = [];

    % Whether to generate a video of RFpooling maps along the horizontal meridian
    options.generateVideoOfRFpoolingMapsAlongHorizontalMeridian (1,1) logical = false;

    % Whether to close previously open figures
    options.closePreviouslyOpenFigures (1,1) logical = true;
end

    % Set flags from key/value pairs
    
    % Mosaic specifiers for selecting a prebaked mRGC mosaic
    rgcMosaicName = options.rgcMosaicName;
    coneMosaicSpecies = options.coneMosaicSpecies;
    opticsSubjectName = options.opticsSubjectName;
    targetVisualSTFdescriptor = options.targetVisualSTFdescriptor;


    %% Close all figures
    close all;

    % Configure a conservative parpool manager. This gives at least 8 GB RAM/core
    % You might need to set this to true if you are hitting a strange
    % parpool crash.  Setting it to true, however, invokes a new parpool
    % open and restore which slows things down in most cases.
    fancyParpoolControl = false;
    if (fancyParpoolControl)
        ASPPManager = AppleSiliconParPoolManager('conservative');
    end

    % Control saving of figures.  We don't want tutorials
    % saving things into the isetbio source tree.
    saveFigures = true;
    figureDir = figExporter.figureDir(mfilename, saveFigures);
    

    %% Load the mRGCMosaic
    % This tutorial only visualizes mosaic structure, so we skip generating
    % the optics. Doing so avoids the Strehl ratio defocus search, which
    % dominates the run time. See t_mRGCMosaicVisualizeWithOptics for
    % visualizations that do require the optics.
    theMRGCMosaic = mRGCMosaic.loadPrebakedMosaic(...
            coneMosaicSpecies, opticsSubjectName, rgcMosaicName, targetVisualSTFdescriptor, ...
            'computeTheMosaicOptics', false);
   

    %% Visualize the mRGCMosaic
    % Lets visualize it together with its input cone mosaic.
    % In this visualization, the gray contours identify the mRGC RF centers, 
    % whereas cones connected to mRGCs are displayed by colored disks.
    theMRGCMosaic.visualize(...
        'identifyInputCones', ~true, ...
        'identifyPooledCones', ~true, ...
        'identifiedConeAperture', 'lightCollectingAreaCharacteristicDiameter', ...
        'plotTitle', 'full mosaic');

    %% Crop the mRGCMosaic to half size
    theMRGCMosaic.cropToSizeAtEccentricity(theMRGCMosaic.sizeDegs*0.5, theMRGCMosaic.eccentricityDegs);

    theMRGCMosaic.visualize(...
        'identifyInputCones', true, ...
        'identifyPooledCones', true, ...
        'identifiedConeAperture', 'lightCollectingAreaCharacteristicDiameter', ...
        'centerSubregionContourSamples', 24, ...
        'plotTitle', 'cropped mosaic');

    %% Visualize the retinal cone pooling for a couple neurons
    % Lets visualize the retinal cone pooling of a couple of neurons.
    % For the first neuron visualization, find a unit with 2 L-cone
    % inputs in its RF center, and located near eccentricityDegs+[-0.5 0]; degs. 
    % 

    targetRGCposition = theMRGCMosaic.eccentricityDegs+ [-0.1 0.1];
    targetCenterConesNum = 2;
    targetCenterConeMajorityType = cMosaic.MCONE_ID;
    theRGCindex(1) = theMRGCMosaic.visualizeRetinalConePoolingRFmapNearPosition(...
       targetRGCposition, targetCenterConesNum, ...
       targetCenterConeMajorityType, ...
       'tickSeparationArcMin', 3);

    targetCenterConesNum = 3;
    targetCenterConeMajorityType = cMosaic.LCONE_ID;
    theRGCindex(2) = theMRGCMosaic.visualizeRetinalConePoolingRFmapNearPosition(...
       targetRGCposition, targetCenterConesNum, ...
       targetCenterConeMajorityType, ...
       'tickSeparationArcMin', 3);

    theMRGCMosaic.visualize(...
        'identifiedConeAperture', 'lightCollectingAreaCharacteristicDiameter', ...
        'identifyInputCones', true, ...
        'identifyPooledCones', true, ...
        'labelRGCsWithIndices', theRGCindex, ...
        'labeledRGCsColor', [0 0 1], ...
        'labeledRGCsLineWidth', 4.0, ...
        'centerSubregionContourSamples', 24, ...
        'backgroundColor', [1 1 1]);

    % Restore previous parpool config
    if (fancyParpoolControl)
        ASPPManager.restoreLastParpoolSize();
    end

end
