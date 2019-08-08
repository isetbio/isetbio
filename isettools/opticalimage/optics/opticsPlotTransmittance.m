function opticsPlotTransmittance(oi)
% Plot spectral transmittance of the optics
%
% Syntax:
%   opticsPlotTransmittance(oi)
%
% Description:
%    Plot the transmittance of the lens and other intervening media.
%
%    This slot is used to store the human macular pigment density. It can
%    also be used for the lens transmittance or the combination of the two.
%
% Inputs:
%    oi - Struct. An optical image structure.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    03/09/18  jnm  Formatting

figNum =  vcSelectFigure('GRAPHWIN');
plotSetUpWindow(figNum);

wave = oiGet(oi, 'wave');
optics = oiGet(oi, 'optics');
transmittance = opticsGet(optics, 'transmittance');
if isempty(transmittance), transmittance = ones(size(wave)); end

plot(wave, transmittance, '-o')

udata.wave = wave; udata.transmittance = transmittance;
set(gca, 'userdata', udata);
xlabel('Wavelength (nm)');
ylabel('Transmittance');
title('Optical transmittance');
grid on

end