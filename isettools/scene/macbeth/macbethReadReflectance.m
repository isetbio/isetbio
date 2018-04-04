function macbethChart = macbethReadReflectance(wave, patchList)
% Read the macbeth surface reflectances into the standard vimage ordering
%
% Syntax:
%	macbethChart = macbethReadReflectance([wave], [patchList]);
%
% Description:
%    The returned variable has the reflectances in the columns but
%    according to the ordering used in vcimage
%
%    The gray series is in macbethChart(:, 4:4:end).
%    The upper left corner (:, 1) is brown, the upper right (:, 21) is cyan
%    The third row is Blue, green, red, yellow, magenta, light blue
% 
%    The code contains examples of function usage. To access, enter 'edit
%    macbethReadReflectance.m' into the Command Window.
%
% Inputs:
%    wave         - (Optional) The wavelength(s). Default is 400:700
%    patchList    - (Optional) The list of patches (colors) to display.
%                   Default is 1:24 (all 24 patch colors).
%
% Outputs:
%    macbethChart - The vimage ordering of macbeth surface reflectances
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    t_SurfaceModels and macbethChartCreate: These are examples of
%    calculating images and other values with the MCC
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    01/30/18  jnm  Formatting

% Examples:
%{
    wave = 400:10:700;
    macbethReflectance = macbethReadReflectance(wave);
    plot(wave, macbethReflectance), xlabel('Wavelength (nm)')
    ylabel('Reflectance');
    grid on
%}
if notDefined('wave'), wave = (400:700); end
if notDefined('patchList'), patchList = 1:24; end

fName = fullfile(isetbioDataPath, 'surfaces', 'macbethChart.mat');
macbethChart = ieReadSpectra(fName, wave);

macbethChart = macbethChart(:, patchList);

end
