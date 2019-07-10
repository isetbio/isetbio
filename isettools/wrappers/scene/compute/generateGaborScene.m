function scene = generateGaborScene(varargin)
% Method to generate an ISETBio scene representing a Gabor stimulus
%
% Syntax:
%   scene = generateGaborScene(varargin])
%
% Description:
%    This function generates an ISETBio scene of a Gabor stimulus based on 
%    the passes stimulus parameters. We use a built-in scene generation 
%    method in ISETBio which generates sinusoidal images. This function 
%    requires a params struct which we generate based on the stimParams info
%
% Optional input key/value pairs:
%   stimParams                 - Stimulus parameters for the Gabor. If not
%                                passed the defaultStimParams is used, whose
%                                fields are defined below.
%
%   presentationDisplay        - Display object, If passed, the returned
%                                scene will be the scene as realized on
%                                that display. In this case the value of 
%                                'pixelsAlongWidthDim' and 'pixelsAlongHeightDim'
%                                will be overwritten based on the passed
%                                display's pixel size
%
% Outputs:
%    scene                     - The ISETBio scene representing the Gabor
%

% History:
%    11/23/18  NPC  ISETBIO TEAM, 2018

%% default stimParams. If a stimParams is passed, it must have these and 
%% only these fields.
defaultStimParams = struct(...
    'spatialFrequencyCyclesPerDeg', 10, ... % 10 cycles/deg
    'orientationDegs', 0, ...               % 45 degrees
    'phaseDegs', 0, ...                     % spatial phase in degrees
    'sizeDegs', 0.5, ...                    % 0.5 x 0.5 degrees
    'sigmaDegs', 0.2/3, ...                 % sigma of Gaussian envelope
    'contrast', 0.6,...                     % 0.6 Michelson contrast
    'meanLuminanceCdPerM2', 40, ...         % mean luminance
    'pixelsAlongWidthDim', 256, ...         % pixels- width dimension
    'pixelsAlongHeightDim', 256 ...         % pixels- height dimension
    );

%% parse input
p = inputParser;
p.addParameter('stimParams', defaultStimParams, @validateStimParamsStruct);
p.addParameter('presentationDisplay', [], @validateDisplayArgument);
p.addParameter('minimumPixelsPerHalfPeriod', [], @isnumeric);
p.parse(varargin{:});

presentationDisplay = p.Results.presentationDisplay;
stimParams = p.Results.stimParams;
minimumPixelsPerHalfPeriod = p.Results.minimumPixelsPerHalfPeriod;

if (~isempty(presentationDisplay))
    [stimParams, presentationDisplay] = updateStimParamsForDisplay(stimParams, presentationDisplay, minimumPixelsPerHalfPeriod);
end

% Transform stimParams into imageHarmonicParams
imageHarmonicParams = struct(...
        'freq', stimParams.sizeDegs * stimParams.spatialFrequencyCyclesPerDeg, ...
        'ang', stimParams.orientationDegs/180*pi, ...
        'ph', stimParams.phaseDegs/180*pi, ...
        'contrast', stimParams.contrast, ...
        'GaborFlag', stimParams.sigmaDegs/stimParams.sizeDegs, ...
        'row', stimParams.pixelsAlongHeightDim , ...
        'col', stimParams.pixelsAlongWidthDim);
    
% Generate a scene using these params
scene = sceneCreate('harmonic', imageHarmonicParams);

% Match viewing distance and mean luminance params
scene = sceneSet(scene, 'wangular', stimParams.sizeDegs);
scene = sceneAdjustLuminance(scene, stimParams.meanLuminanceCdPerM2);

if (~isempty(presentationDisplay))
    % Place the scene at the same distance as the display
    viewingDistanceMeters = displayGet(presentationDisplay, 'viewing distance');
    scene = sceneSet(scene, 'distance', viewingDistanceMeters);
    % Realize the scene into the presentation display
    scene = realizeSceneOnDisplay(scene, presentationDisplay);
end

% Method to validate the passed stimParams struct. The passed stimParams
% struct must have ALL and ONLY the fields of the defaultStimParams
% struct
function vStatus = validateStimParamsStruct(x)

    vStatus = true;
    % Check 1. Input is a struct 
    if (~isstruct(x))
        error('Passed stimParams is not a struct.');
    end
    % Check fields now
    defaultFieldNames = fieldnames(defaultStimParams);
    passedFieldNames = fieldnames(x);
    
    % Check 2. Ensure that we know all the passed field names
    for k = 1:numel(passedFieldNames)
        theFieldName = passedFieldNames{k};
        if (~ismember(theFieldName, defaultFieldNames))
            error('Dont know how to use ''stimParms.%s''.', theFieldName);
        end
    end
    
    % Check 3. Ensure that is all an expected field names have been defined
    for k = 1:numel(defaultFieldNames)
        theFieldName = defaultFieldNames{k};
        if (~ismember(theFieldName, passedFieldNames))
            error('The value of ''stimParms.%s'' has not been defined.', theFieldName);
        end
    end
end

end

% Method to compute the number of pixels for the stimulus given the
% stimulus size, and the viewing distance & pixel size of the display
function [stimParams, presentationDisplay] = updateStimParamsForDisplay(stimParams, presentationDisplay, minimumPixelsPerHalfPeriod)

    if ((~isempty(minimumPixelsPerHalfPeriod)) && (stimParams.spatialFrequencyCyclesPerDeg > 0))
        stimulusHalfPeriodDegrees = 0.5 * 1.0/stimParams.spatialFrequencyCyclesPerDeg;
        displayPixelSizeDegrees = 0.5*displayGet(presentationDisplay, 'deg per dot');
        pixelsPerHalfPeriod = stimulusHalfPeriodDegrees / displayPixelSizeDegrees;

        % We want a minimum of 5 pixels / stimulus half period
        if (pixelsPerHalfPeriod < minimumPixelsPerHalfPeriod)
            % boost the display's resolution to achieve 5 pixels / stimulus half period
            resolutionBoostFactor = minimumPixelsPerHalfPeriod/pixelsPerHalfPeriod;
            oldDPI = displayGet(presentationDisplay, 'dpi');
            newDPI = ceil(oldDPI*resolutionBoostFactor);
            presentationDisplay = displaySet(presentationDisplay, 'dpi', newDPI);
        end
    end
    
    displayPixelSizeDegrees = 0.5*displayGet(presentationDisplay, 'deg per dot');
    % divide by the stimulus size in degrees to get the pixels along the width
    stimParams.pixelsAlongWidthDim = ...
        round(stimParams.sizeDegs/displayPixelSizeDegrees(1));
    stimParams.pixelsAlongHeightDim = ...
        round(stimParams.sizeDegs/displayPixelSizeDegrees(1));
    
end

% Method to validation the passed display argument
function vStatus = validateDisplayArgument(x)

   vStatus = true;
   if isempty(x)
       % OK, we can have an empty argument and we do nothing
   elseif ~isstruct(x)
       % If non-empty it must be a struct object
       error('Input is not an ISETBio display object struct');
   else
       if ~isfield(x, 'type')
           error('Input is not an ISETBio display object struct');
       end
       if ~strcmp(x.type, 'display')
            error('Input is not an ISETBio display object struct');
       end
   end
end