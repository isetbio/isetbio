function fig = CalibratePlotSpectra(calOrCalStruct,fig)
% fig = CalibratePlotSpectra(calOrCalStruct,[fig])
%
% Make a diagnostic plot of the device spectral data and fits in the
% calibration structure.
%
% Can pass figure handle. Returns figure handle.
%
% See also CalibratePlotGamma, CalibratePlotAmbient.
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
rawGammaTable = calStructOBJ.get('rawGammaTable');
S        = calStructOBJ.get('S');
P_device = calStructOBJ.get('P_device');

clf; hold on
if (size(rawGammaTable,2) > 3)
    subplot(1,2,1);
end
hold on
plot(SToWls(S), P_device(:,1),'r');
plot(SToWls(S), P_device(:,2),'g');
plot(SToWls(S), P_device(:,3),'b');
xlabel('Wavelength (nm)', 'Fontweight', 'bold');
ylabel('Power', 'Fontweight', 'bold');
title('Phosphor spectra', 'Fontsize', 13, 'Fontname', 'helvetica', 'Fontweight', 'bold');
axis([380,780,-Inf,Inf]);
if (size(rawGammaTable,2) > 3)
    subplot(1,2,2); hold on
    plot(SToWls(S), P_device(:,4),'r');
    plot(SToWls(S), P_device(:,5),'g');
    plot(SToWls(S), P_device(:,6),'b');
    xlabel('Wavelength (nm)', 'Fontweight', 'bold');
    ylabel('Power', 'Fontweight', 'bold');
    title('Phosphor correction', 'Fontsize', 13, 'Fontname', 'helvetica', 'Fontweight', 'bold');
    axis([380,780,-Inf,Inf]);
end
drawnow;

% calStructOBJ is not needed anymore. So clear it from the memory.
clear 'calStructOBJ'
    