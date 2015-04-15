function displayObject = GenerateIsetbioDisplayObjectFromCalStructObject(displayName, calStructOBJ, varargin)
% displayObject = generateIsetbioDisplayObjectFromCalStructObject(displayName, calStructOBJ, varargin)
%
% Method to generate an isetbio display object with given specifications
%
% 2/20/2015    npc  Wrote skeleton script for xiamao ding
% 2/24/2015    xd   Updated to compute dpi, and set gamma and spds
% 2/26/2015    npc  Updated to employ SPD subsampling 
% 3/1/2015     xd   Updated to take in optional S Vector 
% 3/2/2015     xd   Updated to take in ExtraCalData struct
% 3/9/2015     xd   Updated S Vector behavior
% 4/15/2015    npc  Cleaned up a bit, subsample Svector is now a property of ExtraData
    
    % Input parser to see if optional S vector input exists
    input = inputParser;
    
    % Check that the CalStruct input is indeed a CalStruct
    checkCalStruct = @(x) isa(x, 'CalStruct');

    % Check is ExtraCalData
    checkExtraData = @(x) isa(x, 'ExtraCalData');
    
    addRequired(input, 'displayName', @ischar);
    addRequired(input, 'calStructOBJ', checkCalStruct);
    addRequired(input, 'ExtraData', checkExtraData);
    parse(input, displayName, calStructOBJ, varargin{:});
    
    % Assemble filename for generated display object
    displayFileName = sprintf('%s.mat', displayName);
    
    % Generate a display object
    displayObject = displayCreate;

    % Set the display's name to the input parameter displayName
    displayObject = displaySet(displayObject, 'name', displayFileName);

    % Get the wavelength sampling and the SPD from the CalStructOBJ 
    S = calStructOBJ.get('S');
    spd = calStructOBJ.get('P_device');

    if (~isempty(input.Results.ExtraData.subSamplingSvector))
        % Validate that the subSamplingSvector is within range of the original S vector
        validateSVector(calStructOBJ.get('S'), input.Results.ExtraData.subSamplingSvector);
        % SubSample the SPDs 
        newS = input.Results.ExtraData.subSamplingSvector;              
        lowPassSigmaInNanometers = 4;        
        maintainTotalEnergy = true;
        showFig = false;
        if (exist('subSampleSPDs', 'file'))
            fprintf('Will subsample SPDs with a resolution of %d nm\n', newS(2));
            [subSampledWave, subSampledSPDs] = subSampleSPDs(S, spd, newS, lowPassSigmaInNanometers, maintainTotalEnergy, showFig);
        else
            error('Function ''subSampleSPDs'', is a function of the BLIlluminationDiscrimination project, but this project is not on your path');
        end
        % Set the display object's SPD to the subsampled versions
        displayObject = displaySet(displayObject, 'wave', subSampledWave);
        displayObject = displaySet(displayObject, 'spd', subSampledSPDs);
    else
        fprintf('Will not subsample SPDs\n');
        % Set the display object's SPD to the original versions
        displayObject = displaySet(displayObject, 'wave', SToWls(S));
        displayObject = displaySet(displayObject, 'spd', spd);
    end

    % Get the display's gamma table
    displayObject = displaySet(displayObject, 'gTable', calStructOBJ.get('gammaTable'));

    % Get the display resolution in dots (pixels) per inch
    m = calStructOBJ.get('screenSizeMM');
    p = calStructOBJ.get('screenSizePixel');

    m = m/25.4;
    mdiag = sqrt(m(1)^2 + m(2)^2);
    pdiag = sqrt(p(1)^2 + p(2)^2);
    dpi = pdiag / mdiag;

    displayObject = displaySet(displayObject, 'dpi', dpi);

    % Use the viewing distance obtained from the ExtraData Struct
    dist = input.Results.ExtraData.distance;
    displayObject = displaySet(displayObject, 'viewing distance', dist);
    
    % Save display object to file
    fprintf('Saving new display object (''%s'').\n', displayName);
    d = displayObject;
    save(displayFileName, 'd');   
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