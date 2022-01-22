% Generate a PTB calibration structure from an ISETBIO display object.
% 
% Synopsis:
%    PTBcal = ptb.GeneratePTCalStructFromIsetbioDisplayObject(display)
%
% Description:
%    Produce a PTB calibration struct from an ISETBio display object.
%
%    Note that the PTB power units convention is power per wavelength band, 
%    while ISETBio is power per nm. This routine does the conversion as it extracts
%    spectra from ISETBio and puts them into the PTB calibration structure.
%    The inverse routine, ptb.GenerateIsetbioDisplayObjectFromPTBCalStruct,
%    does the same conversion in the other direction. 
%
%    This routine assumes that the display object is describing a 3 primary
%    device. If there are additional primary spectra in the display object, they
%    are ignored.
%
%    Note unfortunate choice of PT rather than PTB in the name of this
%    function. But changing it now will likely break extant code and
%    doesn't seem worth it.
%
% Inputs:
%    display           - The ISETBio display object.
%
% Outputs:
%    PTBCal            - The PTB calibration structure.
%
% Optional key/value pairs:
%   None.
%
% See also: ptb.GenerateIsetbioDisplayObjectFromPTBCalStruct,
%           ptb.GeneratePsychToolboxCalStruct, ptb.GenerateEmptyPTBStruct,
%           generateCustomDisplay.
%

% History:
%   6/25/15  dhb  Modify to get display size out of ISETBIO object and use it.
%   1/16/22  dhb  Do the ISETBio -> PTB power unit convention conversion.
%                 And add some comments.

function PTBcal = GeneratePTCalStructFromIsetbioDisplayObject(display)

    % Start with a totally empty PTBcal
    PTBcal = ptb.GenerateEmptyPTBcalStruct();
    
    % Update key properties
    PTBcal = updateDisplayDescription(PTBcal, display);
    PTBcal = updateSpectralParams(PTBcal, display);
    PTBcal = updateGammaParams(PTBcal, display);  
end
    
function PTBcal = updateGammaParams(oldPTBcal, display)
    [gammaTable, gammaInput] = retrieveGammaTable(display);
    
    PTBcal = oldPTBcal;
    PTBcal.gammaInput = gammaInput';
    PTBcal.gammaTable = gammaTable;
    PTBcal.nDevices = size(gammaTable,2);
end

function PTBcal = updateSpectralParams(oldPTBcal, display)
    
    [wave, spd] = retrievePrimaries(display);
    spectralSamples  = size(wave,1);
    
    PTBcal = oldPTBcal;
    PTBcal.describe.S   = WlsToS(wave);
    PTBcal.S_ambient    = PTBcal.describe.S;
    PTBcal.P_device     = spd*PTBcal.describe.S(2);
    PTBcal.P_ambient    = displayGet(display,'ambient spd')*PTBcal.describe.S(2);
    PTBcal.T_ambient    = eye(spectralSamples);
    PTBcal.T_device     = eye(spectralSamples); 

    % For an ISETBio display, this is always 1.
    PTBcal.nPrimaryBases = 1;
end

function PTBcal = updateDisplayDescription(oldPTBcal,display)   
    dotsPerMeter = displayGet(display, 'dots per meter');
    screenSizeMM = 1000.0*displayGet(display,'size');
    screenSizePixels = round(screenSizeMM/1000*dotsPerMeter);
    
    PTBcal = oldPTBcal;
    PTBcal.describe.displayDescription.screenSizePixel = screenSizePixels;
    PTBcal.describe.displayDescription.screenSizeMM = screenSizeMM;
end


function [gammaTable, gammaInput] = retrieveGammaTable(display)
    % Gamma table, remove 4-th primary, if it exists
    gammaTable = displayGet(display, 'gTable');
    if (size(gammaTable,2) > 3)
        gammaTable = gammaTable(:,1:3);
    end
    gammaInput = linspace(0,1,size(gammaTable,1));
end

function [wave, spd] = retrievePrimaries(display)
    % Remove 4-th primary, if it exists, for testing purposes.
    wave = displayGet(display, 'wave');
    spd  = displayGet(display, 'spd primaries');
    if (size(spd ,2) > 3)
        spd = spd(:,1:3);
    end
end


