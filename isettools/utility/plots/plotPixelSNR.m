function [uData, figHdl] = plotPixelSNR(sensor)
% Graph the pixel SNR over the pixel response range
%
% Syntax:
%   [uData, figHdl] = plotPixelSNR(sensor)
%
% Description:
%    Three curves are generated. They show the following information:
%       1) The total pixel SNR.
%       2) The SNR limits from Shot Noise
%       3) The SNR limits from Read Noise
%
%    This routine usese the currently selected ISA structure to retrieve
%    all of the data and properties used for plotting. Perhaps it should
%    take an ISA argument.
%
% Inputs:
%    sensor - The sensor object
%
% Outputs:
%    uData  - The user data structure
%    figHdl - The figure handle
%
% Notes:
%    * [Note: JNM - The function pixelSNR does not exist!]
%    * [Note: JNM - legend using a number to specify the location is not
%      supported. What do we need to change about this?]
%

% History:
%    xx/xx/12       (c) Imageval Consulting, LLC, 2012
%    12/11/17  jnm  Formatting
%

% Example:
%{
    plotPixelSNR;
    [uData, g] = plotPixelSNR(vcGetObject('sensor'));
%}

if notDefined('sensor'), sensor = vcGetObject('sensor'); end
% SNRread can come back as Inf if there is no read noise.
[SNR, volts, SNRshot, SNRread] = pixelSNR(sensor);

uData.volts = volts;
uData.snr = SNR;
uData.snrShot = SNRshot;
uData.snrRead = SNRread;

figHdl = vcNewGraphWin;

p = semilogx(volts, SNR, 'k-');
hold on;
set(p, 'linewidth', 2);
p = semilogx(volts, SNRshot, 'r-', volts, SNRread, 'g-');
set(p, 'linewidth', 1);
hold off
grid on; 
xlabel('Signal (V)'); 
ylabel('SNR (db)');
title('Pixel SNR over response range');
legend('Total pixel SNR', 'Shot noise SNR', 'Read noise SNR', 2);

% Attach data to the figure
set(figHdl, 'Userdata', uData);
return;
