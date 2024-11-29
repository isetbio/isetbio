% RefitCalLinMod
%
% Refit the calibration linear model.
%
% 3/27/02  dhb            Wrote it.
% 9/26/08  dhb, ijk, tyl  Simplify naming possibilities. 
% 9/27/08  dhb            Clearer defaults for prompts.  Pass number of levels to dacsize routine.
% 2/15/10  dhb            Plot all components of gamma functions.
% 5/28/10  dhb            Add yoked fitting routine to calls.  Should have no effect when yoked isn't set, but 
%                         do the right thing when it is.
% 6/5/10   dhb            Make it work when rawGammaInput has multiple columns.  Use plot subroutines.
% 5/28/14  npc            Modifications for accessing calibration data using a @CalStruct object.

% Enter load code
defaultFileName = 'StereoLCDLeft';
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
oldN            = calStructOBJ.get('nPrimaryBases');

% Provide information about gamma measurements
% This is probably not method-independent.
fprintf('Gamma measurements were made at %g levels\n',...
	size(rawGammaInput,1));
fprintf('Gamma table available at %g levels\n',...
	size(gammaInput,1));

% Get new fit type
fprintf('Old linear model fit was with %g components\n',oldN);
newN = input(sprintf('Enter new number of components: [%d]: ',oldN));
if (~isempty(newN))
    calStructOBJ.set('nPrimaryBases', newN)
end


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





