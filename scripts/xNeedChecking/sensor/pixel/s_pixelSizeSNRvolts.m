%% s_pixelSizeSNRvolts
%
%  UNDER DEVELOPMENT
%  How pixel properties influence signal-to-noise ratio (SNR)
%
% The pixel SNR is specified as
%
%     10 log10 ( sigPower / noisePower )
%
% with definitions
%
%  sigPower:    The signal power at any voltage is the square of the number
%  of electrons at that voltage.  Signal variance is Poisson and thus the
%  signal variance equals the number of electrons at each voltage (i.e.,
%  the square root of the signal power).
%
%  noisePower:  The noise power is the sum of the signal variance and the
%  read noise variance. The read noise variance is a parameter of the technology.
%
% For each technology, the SNR varies as a function of the pixel voltage
% level.  The precise SNR depends on factors like the dark voltage. In
% general, as the pixel voltage  increases, the SNR also increases. 
%
% To compute pixel SNR, we use pixel signal and noise specified in units of
% electrons (see pixelSNR).  This is important because the Poisson
% character of the noise is only true in units of electrons, but not in
% units of volts (the Poisson distribution is not invariant with respect to
% scaling).
%
%
% See also: pixelSizeSNRluxsec
%
% Copyright ImagEval Consultants, LLC, 2005.

% You should set the parameters here according to your technology
% properties

integrationTime = 0.010;                      % Sec
pixelSize    = [2 4 6 9 10]*1e-6;             % Pixel size in meters
readNoiseSD  = [5 4 3 2 1]*1e-3;              % std dev in volts
voltageSwing = [.7 1.2 1.5 2 3];              % voltage swing
darkVoltage  = [1 1 1 1 1]*1e-3;              % Volts per sec

% Now you can run this code.
sensor = sensorCreate('monochrome');                %Initialize
pixel = sensorGet(sensor,'pixel');
sensor = sensorSet(sensor,'integrationTime',integrationTime);

for ii=1:length(pixelSize)
    pixel = pixelSet(pixel,'size',[pixelSize(ii),pixelSize(ii)]);
    pixel = pixelSet(pixel,'readNoiseSTDvolts',readNoiseSD(ii));
    pixel = pixelSet(pixel,'voltageSwing',voltageSwing(ii));
    pixel = pixelSet(pixel,'darkVoltage',darkVoltage(ii));

    sensor = sensorSet(sensor,'pixel',pixel);
    [a,b] =  pixelSNR(sensor);
    SNR{ii} = a; volts{ii} = b;
end

% The data were saved in the cell arrays SNR{} and volts{}.  Here is a
% summary plot. 
figure;
txt = sprintf('Pixel size\n');
c = {'r','g','b','c','m','y','k'};
for ii=1:length(SNR)
    semilogx(volts{ii},SNR{ii},['-',c{ii}])
    hold on;
    newText = sprintf('%.0f um\n',pixelSize(ii)*10^6);
    txt = addText(txt,newText);
end
plotTextString(txt,'ur');
hold off
xlabel('Pixel voltage'), ylabel('SNR (db)');
title('SNR vs. voltage')
grid on

% End of Script