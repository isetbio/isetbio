function displayObject = GenerateIsetbioDisplayObjectFromPTBCalStruct(displayName, calStruct, ExtraCalData, saveDisplayObject)
% Generate an isetbio display object with given specifications.
% 
% Synopsis:
%     displayObject = ptb.GenerateIsetbioDisplayObjectFromCalStructObject(displayName, calStruct, ExtraCalData, saveDisplayObject)
%
% Description:
%    Convert a PTB calibration structure to an ISETBio display object, with
%    optional wavelength resampling specified in a property of input
%    argument ExtraCalData.
%
%    If ExtraCalData has a subsampling S vector in it, the spectra are
%    subsampled with a lowpass filter consisting of a Gaussian kernal of 4
%    nm sd.  That may be oversmoothing for some choices. An enterprising
%    person could update this routine to make the sd an optional key/value
%    pair for this routine, and then update this comment as well as the
%    optional key/value pairs section below.
%
%    Note that the PTB power units convention is power per wavelength band,
%    while ISETBio is power per nm. This routine does the conversion as it
%    extracts spectra from the PTB calibration structure and puts them into
%    the isetbio display object. The inverse routine,
%    ptb.GeneratePTCalStructFromIsetbioDisplayObject, does the same
%    conversion in the other direction.
%
% Inputs:
%   displayName        - Name to give the created display
%   calStruct          - PTB calibration structure to convert.
%   ExtraCalData       - Object whose properties contain extra data not in
%                        the PTB calibration structure but needed to
%                        produce a proper ISETBio display object.  See
%                        ptb.ExtraCalData.
%   saveDisplayObject  - The display obect is saved in displayName.mat if
%                        this is true, and not if it is false.
%
% Outputs:
%   displayObject      - ISETBio display object.
%
% Optional key/value pairs:
%   None. But it might be clever to add a 'lowPassSigmaInNanometers' key
%   value pair as the default of 4 nm won't always be the best choice. An
%   alternative approach would be to add a property to the ExtraCalData
%   class and use that.
%
% See also: ptb.ExtraCalData, ptb.SubSampleSPDs, 
%           ptb.GeneratePTCalStructFromIsetbioDisplayObject,
%           ptb.GeneratePsychToolboxCalStruct, ptb.GenerateEmptyPTBStruct,
%           generateCustomDisplay.

% History:
%   2/20/2015    npc  Wrote skeleton script for xiamao ding
%   2/24/2015    xd   Updated to compute dpi, and set gamma and spds
%   2/26/2015    npc  Updated to employ SPD subsampling 
%   3/1/2015     xd   Updated to take in optional S Vector 
%   3/2/2015     xd   Updated to take in ExtraCalData struct
%   3/9/2015     xd   Updated S Vector behavior
%   4/15/2015    npc  Cleaned up a bit, subsample Svector is now a property of ExtraCalData
%   4/15/2015    npc  Added input arg, to control whether to save the generated isetbio display object  
%   6/25/15      dhb  Set ISETBIO display size field

    % Check is ExtraCalData
    checkExtraCalData = @(x) isa(x, 'ptb.ExtraCalData');
    
    % Input parser to check validity of inputs
    input = inputParser;
    addRequired(input, 'displayName', @ischar);
    addRequired(input, 'calStruct', @isstruct);
    addRequired(input, 'ExtraCalData', checkExtraCalData);
    addRequired(input, 'saveDisplayObject', @islogical);
    parse(input, displayName, calStruct, ExtraCalData, saveDisplayObject);

    % Set this to true for some diagnostic output
    showFig = false;
    
    % Assemble filename for generated display object
    displayFileName = sprintf('%s.mat', displayName);
    
    % Generate a display object
    displayObject = displayCreate;

    % Set the display's name to the input parameter displayName
    displayObject = displaySet(displayObject, 'name', displayFileName);

    % Get the wavelength sampling and channel spds, and ambient spd from
    % the input calStruct. The spds are in PTB units of power per
    % wavelength band.
    S = calStruct.describe.S;
    spd = calStruct.P_device;
    ambient = calStruct.P_ambient;
    
    if (~isempty(input.Results.ExtraCalData.subSamplingSvector))
        % Subsample and set the spectra
        if (showFig)
            fprintf('Will subsample SPDs with a resolution of %d nm\n', input.Results.ExtraCalData.subSamplingSvector(2));
        end

        % Validate that the subSamplingSvector is within range of the original S vector
        validateSVector(S, input.Results.ExtraCalData.subSamplingSvector);
        
        % SubSample the SPDs 
        newS = input.Results.ExtraCalData.subSamplingSvector;              
        lowPassSigmaInNanometers = 4;        

        % The ptb.SubSampleSPDs routine converts the power units.
        [subSampledWave, subSampledSPDs] = ptb.SubSampleSPDs(S, spd, newS, lowPassSigmaInNanometers, showFig);
        [~, subSampledAmbient] = ptb.SubSampleSPDs(S, ambient, newS, lowPassSigmaInNanometers,  showFig);

        % Set the display object's SPD to the subsampled versions
        displayObject = displaySet(displayObject, 'wave', subSampledWave);
        displayObject = displaySet(displayObject, 'spd', subSampledSPDs);
        displayObject = displaySet(displayObject, 'ambient spd', subSampledAmbient);
    else
        % Set the display object's SPD to the original versions.  Need to
        % do the unit conversion right here, right now. In the parallel
        % sets just above we don't, because the conversion is done inside
        % of routine ptbSubSampleSPDs.
        if (showFig)
            fprintf('Will not subsample SPDs\n');
        end
        displayObject = displaySet(displayObject, 'wave', SToWls(S));
        displayObject = displaySet(displayObject, 'spd', spd/S(2));
        displayObject = displaySet(displayObject, 'ambient spd', ambient/S(2));
    end

    % Get the display's gamma table.
    gammaTable = calStruct.gammaTable;
    gammaLength = size(gammaTable,1);
    displayObject = displaySet(displayObject, 'gTable', gammaTable);

    % Get the display resolution in dots (pixels) per inch
    m = calStruct.describe.displayDescription.screenSizeMM;
    p = calStruct.describe.displayDescription.screenSizePixel;
    m = m/25.4;
    mdiag = sqrt(m(1)^2 + m(2)^2);
    pdiag = sqrt(p(1)^2 + p(2)^2);
    dpi = pdiag / mdiag;
    displayObject = displaySet(displayObject, 'dpi', dpi);
    
    % Set the display size 
    displayObject = displaySet(displayObject,'size',calStruct.describe.displayDescription.screenSizeMM/1000);

    % Use the viewing distance obtained from the ExtraCalData Struct
    dist = input.Results.ExtraCalData.distance;
    displayObject = displaySet(displayObject, 'viewing distance', dist);
    
    if (input.Results.saveDisplayObject)
        % Save display object to file
        fprintf('Saving new display object (''%s'').\n', displayName);
        d = displayObject;
        save(displayFileName, 'd');   
    end
    
end

function validateSVector(oldS, newS)  
    % Check that newS fits S vector parameters
    SVecAttribute = {'size', [1,3]};
    SVecClass = {'double'};
    validateattributes(newS, SVecClass, SVecAttribute)
    
    % Check that newS is within range of oldS
    newWave = SToWls(newS);
    oldWave = SToWls(oldS);
    if newS(1) < oldS(1)
        error('S Vector starts at lower nm than original');
    elseif newWave(end) > oldWave(end)
        error('S Vector ends at higher nm than original');
    end
end