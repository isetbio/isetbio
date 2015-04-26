function varargout = v_PTBcalStructToIsetbioDisplayObjectAndBack(varargin)
%
% Validate conversion of PTBcalStruct to isetbio display object and back.
%
    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)

    %% Initialize ISET
    close all; ieInit;
    
    %% Some informative text
    UnitTest.validationRecord('SIMPLE_MESSAGE', 'Conversion of PTB cal struct to isetbio display object - This is still in progress.');
    
    % Load PTB calibration data for the SONY PVM2541 display
    calDir = fullfile(isetbioRootPath, 'isettools', 'data', 'ptbcal');
    calFile = 'SonyPVM2541';
    PTBcal  = LoadCalFile(calFile, [], calDir);
    

    % Generate isetbio display object
    extraData = ptb.ExtraCalData;
    extraData.distance = 0.764;

    % Generate an isetbio display object to model the display used to obtain the calibration data
    saveDisplayObject = false;
    isetbioDisplayObject = ptb.GenerateIsetbioDisplayObjectFromPTBCalStruct('SonyPVM2541', PTBcal, extraData, saveDisplayObject);
    reconstructedPTBcal  = ptb.GeneratePTCalStructFromIsetbioDisplayObject(isetbioDisplayObject);
    
    % Start comparisons
    % Display resolution (along horizontal dimension)
    isetbioDotsPerMeter       = displayGet(isetbioDisplayObject, 'dots per meter')
    originalPTBdotsPerMeterX   = PTBcal.describe.displayDescription.screenSizePixel(1)/(PTBcal.describe.displayDescription.screenSizeMM(1)/1000)
    reconstructedPTBdotsPerMeterX = reconstructedPTBcal.describe.displayDescription.screenSizePixel(1)/(reconstructedPTBcal.describe.displayDescription.screenSizeMM(1)/1000)
    
    isetbioDotsPerMeter       = displayGet(isetbioDisplayObject, 'dots per meter')
    originalPTBdotsPerMeterY   = PTBcal.describe.displayDescription.screenSizePixel(2)/(PTBcal.describe.displayDescription.screenSizeMM(2)/1000)
    reconstructedPTBdotsPerMeterY= reconstructedPTBcal.describe.displayDescription.screenSizePixel(2)/(reconstructedPTBcal.describe.displayDescription.screenSizeMM(2)/1000)
    
    
    %% Plot
    if (runTimeParams.generatePlots)
    end
    
end