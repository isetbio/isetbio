% RefitCalGamma
%
% Refit the gamma data from a calibration.
%
% See also CalDemo, CalibrateFitGamma, SetGammaMethod, GamutToSettings
%
% 3/27/02  dhb      Wrote it.
% 8/26/03  dhb, jl  Allow changing dacsize of calibraiton.  Useful for dropping from 10 to 8.
% 2/2/05   dhb, dam Ask for filename to save to, rather than automatically overwrite.
% 9/26/08  dhb, ijk, tyl  Simplify naming possibilities.  Add some better help.
% 9/27/08  dhb      Fix up dacsize fitting.
%                   Clearer prompts (show default values).
% 11/19/09 dhb      Added crtSumPow option.  This works great for our FrontRoom monitor, which
%                   was not well fit by the traditional gamma model.  The work is done in
%                   function CalibrateFitGamma.  See comments there.
% 11/20/09 dhb      More terms in crtSumPow.
% 3/07/10  dhb      Call CalibrateFitLinMod as well.
% 3/08/10  dhb      Update list of fit type options.
% 5/28/10  dhb      Add yoked fitting routine to calls.  Should have no effect when yoked isn't set, but do the right thing when it is.
% 6/5/10   dhb      Update type list provided to user.
%          dhb      Better plots, using plot subroutines.
% 5/26/12  dhb      Added ability to use raw measured data as the fit gamma table.  See comment where that's done below.
% 5/28/14  npc      Modifications for accessing calibration data using a @CalStruct object.

% Enter load code
defaultFileName = 'PTB3TestCal';
thePrompt = sprintf('Enter calibration filename [%s]: ',defaultFileName);
newFileName = input(thePrompt,'s');
if (isempty(newFileName))
    newFileName = defaultFileName;
end
fprintf(1,'\nLoading from %s.mat\n',newFileName);
cal = LoadCalFile(newFileName);
fprintf('Calibration file %s read\n\n',newFileName);

% Specify @CalStruct object that will handle all access to the calibration data.
[calStructOBJ, inputArgIsACalStructOBJ] = ObjectToHandleCalOrCalStruct(cal);
% Clear cal, so fields are accessed only via get and set methods of calStruct.
clear 'cal'

% Print out some information from the calibration.
DescribeMonCal(calStructOBJ);

% Get necessary calibration data
rawGammaInput   = calStructOBJ.get('rawGammaInput');
rawGammaTable   = calStructOBJ.get('rawGammaTable');
gammaInput      = calStructOBJ.get('gammaInput');
gammaTable      = calStructOBJ.get('gammaTable');
oldType         = calStructOBJ.get('gamma.fitType');
oldDacsize      = calStructOBJ.get('dacsize');

% Provide information about gamma measurements
% This is probably not method-independent.
fprintf('Gamma measurements were made at %g levels\n',...
	size(rawGammaInput,1));
fprintf('Gamma table available at %g levels\n',...
	size(gammaInput,1));

% Get new fit type
fprintf('Old gamma fit type was: %s\n',oldType);
fprintf('Possible fit types are defined by routine CalibrateFitGamma\n');
fprintf('See "help CalibrateFitGamma for most up to date options (except for rawdata)\n');
fprintf('Current (May 2012) options are:\n');
fprintf('\tsimplePower\n');
fprintf('\tcrtLinear\n');
fprintf('\tcrtPolyLinear\n');
fprintf('\tcrtGamma\n');
fprintf('\tcrtSumPow\n');
fprintf('\tbetacdf\n');
fprintf('\tsigmoid\n');
fprintf('\tweibull\n');
fprintf('\trawdata\n');

newType = GetWithDefault('Enter new fit type',oldType);
calStructOBJ.set('gamma.fitType', newType);

if (strcmp(newType,'betacdf'))
    useweight = GetWithDefault('Enter weighting constant (-1 -> no weighting)',0.01);
    calStructOBJ.set('gamma.useweight', useweight);
end

% Get new number of bits
fprintf('Old DAC bits was: %d\n',oldDacsize);
newDacsize = input(sprintf('Enter new dacsize: [%d]: ',oldDacsize));
if (~isempty(newDacsize))
    calStructOBJ.set('dacsize', newDacsize);
end

% Now refit
CalibrateFitLinMod(calStructOBJ);
CalibrateFitYoked(calStructOBJ);

% Switch about whether to call standard routine
fitType = calStructOBJ.get('gamma.fitType');
switch (fitType)
    % Literally use the measured data as the fit data.
    % You only want to do this under very unusual circumstances.
    % The reason I'm putting it in is to deal with a DLP projector
    % whose filter wheel I ripped out, and for which the light output
    % is highly non-monotonic with input.  So, I'll measure at every
    % possible input settings and then use exhaustive search of the
    % gamma table to invert.  This is where the output of the standard
    % calibration program pushes the raw data into the typical field.
    % 
    % This is probably a fairly fragile bit of code and should only be used
    % with caution and knowledge aforethought.
    %
    % I put this here rather than in CalibrateFitGamma to avoid a lot of
    % massaging done by that routine, which we do not want.
    case 'rawdata'
        rawGammaInput   = calStructOBJ.get('rawGammaInput');
        rawGammaTable   = calStructOBJ.get('rawGammaTable');
        nDevices        = calStructOBJ.get('nDevices');
        nPrimaryBases   = calStructOBJ.get('nPrimaryBases');
        if (size(rawGammaTable,2) ~= nDevices*nPrimaryBases)
            error('Dimensions of raw data are not correct, given number of devices and linear model size')
        end
        gammaInput = rawGammaInput;
        gammaTable = rawGammaTable;
        gammaTable(gammaTable < 0) = 0;
        gammaTable(gammaTable > 1) = 1;
        calStructOBJ.set('gammaInput', gammaInput);
        calStructOBJ.set('gammaTable', gammaTable);

    % Fit the measured data using standard PTB methods
    otherwise 
        dacsize = calStructOBJ.get('dacsize');
        CalibrateFitGamma(calStructOBJ,2^dacsize);
end

% Put up a plot of the essential data
CalibratePlotSpectra(calStructOBJ,figure(1));
CalibratePlotGamma(calStructOBJ,figure(2));
drawnow;

% Option to save the refit file
saveIt = input('Save new fit data (0->no, 1->yes)? [0]: ');
if (isempty(saveIt))
	saveIt = 0;
end
if (saveIt)
    % Prompt for new file name if we're saving to a name.
    defaultFileName = newFileName;
    thePrompt = sprintf('Enter calibration filename [%s]: ',defaultFileName);
    saveFileName = input(thePrompt,'s');
    if (isempty(saveFileName))
        saveFileName = defaultFileName;
    end
    fprintf(1,'\nSaving to %s.mat\n',saveFileName);
    SaveCalFile(calStructOBJ.cal,saveFileName);
end



