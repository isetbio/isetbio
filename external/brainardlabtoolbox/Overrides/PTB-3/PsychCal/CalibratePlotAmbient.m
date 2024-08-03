function fig = CalibratePlotAmbient(calOrCalStruct,fig)
% fig = CalibratePlotAmbient(calOrCalStruct,[fig])
%
% Make a diagnostic plot of the ambient spectral data
%
% Can pass figure handle. Returns figure handle.
%
% See also CalibratePlotGamma, CalibratePlotSpectra.
%
% 6/5/10   dhb  Wrote it.
% 5/08/14  npc  Modifications for accessing calibration data using a @CalStruct object.
%               The first input argument can be either a @CalStruct object (new style), or a cal structure (old style).
%               Passing a @CalStruct object is the preferred way because it results in 
%               (a) less overhead (@CalStruct objects are passed by reference, not by value), and
%               (b) better control over how the calibration data are accessed.


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
S  = calStructOBJ.get('S');
P_ambient = calStructOBJ.get('P_ambient');

clf; hold on
plot(SToWls(S), P_ambient(:,1),'k');
xlabel('Wavelength (nm)', 'Fontweight', 'bold');
ylabel('Power', 'Fontweight', 'bold');
title('Ambient spectra', 'Fontsize', 13, 'Fontname', 'helvetica', 'Fontweight', 'bold');
axis([380,780,-Inf,Inf]);
drawnow;
