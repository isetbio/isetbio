% RefitCalYoked
%
% Refit a calibration file to make gamma curves based on the raw yoked measurements.
% This requires, of course, that the yoked measurements were taken during calibration,
% which in turn requires the correct settings.
%
% 4/29/10  dhb, kmo, ar  Wrote it.
% 5/25/10  dhb, ar       Yoked field in describe now set elsewhere.
% 5/28/10  dhb           Update to swing on field yokedmethod.
% 6/5/10   dhb           Use plot subroutines.
% 5/28/14  npc           Modifications for accessing calibration data using a @CalStruct object.

% Enter load code
defaultFileName = 'HDRFrontMaterial';
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
gammaInput      = calStructOBJ.get('gammaInput');
oldMethod       = calStructOBJ.get('yokedmethod');

% Provide information about gamma measurements
% This is probably not method-independent.
fprintf('Gamma measurements were made at %g levels\n',...
	size(rawGammaInput,1));
fprintf('Gamma table available at %g levels\n',...
	size(gammaInput,1));

% Get yoked method
newYokedMethod = input(sprintf('Enter yoked method (0 for not yoked): [%d]: ',oldMethod));
if (isempty(newYokedMethod))
	calStructOBJ.set('yokedmethod', newYokedMethod);
end

% Fit
CalibrateFitLinMod(calStructOBJ);
CalibrateFitYoked(calStructOBJ);
dacsize = calStructOBJ.get('dacsize');
CalibrateFitGamma(calStructOBJ,2^dacsize);

% Put up a plot of the essential data
CalibratePlotSpectra(calStructOBJ,figure(1));
CalibratePlotGamma(calStructOBJ,figure(2));

% Option to save the refit file
saveIt = input('Save new fit data (0->no, 1 -> yes)? [0]: ');
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