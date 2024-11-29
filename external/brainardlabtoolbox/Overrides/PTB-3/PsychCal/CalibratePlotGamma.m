function fig = CalibratePlotGamma(calOrCalStruct,fig)
% fig = CalibratePlotGamma(calOrCalStruct,[fig])
%
% Make a diagnostic plot of the gamma data and fits in the
% calibration structure.
%
% Can pass figure handle. Returns figure handle.
%
% See also CalibratePlotSpectra, CalibratePlotAmbient.
%
% 6/5/10  dhb  Wrote it.
% 5/08/14 npc  Modifications for accessing calibration data using a @CalStruct object.
%              The first input argument can be either a @CalStruct object (new style), or a cal structure (old style).
%              Passing a @CalStruct object is the preferred way because it results in 
%              (a) less overhead (@CalStruct objects are passed by reference, not by value), and
%              (b) better control over how the calibration data are accessed.


% Optional figure open
if (nargin < 2 || isempty(fig))
    fig = figure;
end

% Specify @CalStruct object that will handle all access to the calibration data.
[calStructOBJ, inputArgIsACalStructOBJ] = ObjectToHandleCalOrCalStruct(calOrCalStruct);
if (~inputArgIsACalStructOBJ)
     % The input (calOrCalStruct) is a cal struct. Clear it to avoid  confusion.
    clear 'calOrCalStruct';
end
% From this point onward, all access to the calibration data is accomplised via the calStructOBJ.

% Extract needed data
rawGammaInput = calStructOBJ.get('rawGammaInput');
rawGammaTable = calStructOBJ.get('rawGammaTable');
gammaInput    = calStructOBJ.get('gammaInput');
gammaTable    = calStructOBJ.get('gammaTable');

clf;
if (size(rawGammaTable,2) > 3)
    subplot(1,2,1);
end
hold on
if (size(rawGammaInput,2) == 1)
    plot(rawGammaInput, rawGammaTable(:,1),'r+');
    plot(rawGammaInput, rawGammaTable(:,2),'g+');
    plot(rawGammaInput, rawGammaTable(:,3),'b+');
else
    plot(gammaInput(:,1), gammaTable(:,1),'r+');
    plot(gammaInput(:,2), gammaTable(:,2),'g+');
    plot(gammaInput(:,3), gammaTable(:,3),'b+');
end
xlabel('Input value', 'Fontweight', 'bold');
ylabel('Normalized output', 'Fontweight', 'bold');
title('Gamma functions', 'Fontsize', 13, 'Fontname', 'helvetica', 'Fontweight', 'bold');
plot(gammaInput, gammaTable(:,1),'r');
plot(gammaInput, gammaTable(:,2),'g');
plot(gammaInput, gammaTable(:,3),'b');
axis([0 1 0 1.2]);
if (size(rawGammaTable,2) > 3)
    subplot(1,2,2); hold on
    if (size(rawGammaInput,2) == 1)
        plot(rawGammaInput, rawGammaTable(:,4),'r+');
        plot(rawGammaInput, rawGammaTable(:,5),'g+');
        plot(rawGammaInput, rawGammaTable(:,6),'b+');
    else
        plot(rawGammaInput(:,1), rawGammaTable(:,4),'r+');
        plot(rawGammaInput(:,2), rawGammaTable(:,5),'g+');
        plot(rawGammaInput(:,3), rawGammaTable(:,6),'b+');
    end
    xlabel('Input value', 'Fontweight', 'bold');
    ylabel('Normalized output', 'Fontweight', 'bold');
    title('Gamma correction', 'Fontsize', 13, 'Fontname', 'helvetica', 'Fontweight', 'bold');
    hold on
    plot(gammaInput, gammaTable(:,4),'r');
    plot(gammaInput, gammaTable(:,5),'g');
    plot(gammaInput, gammaTable(:,6),'b');
    axis([0 1 -1.2 1.2]);
end
drawnow;

