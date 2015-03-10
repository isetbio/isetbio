function voltImages = sensorComputeSamples(sensorNF,nSamp,noiseFlag,showBar)
%Computing multiple noise samples of the sensor voltage image
%
%  voltImages = sensorComputeSamples(sensorNF,[nSamp = 10],[noiseFlag=2],[showBar = 1])
%
% Compute multiple noisy samples of the sensor voltage image.  
%
% sensorNF:  Sensor struct containing noise free voltages 
% nSamp:     Number of samples, default = 10
% noiseFlag: photon noise only (1)  all noise (2) Default = 2
% showBar:   Show waitbar (1) or not (0).  Default from ieSession
%
% The voltImages is a 3d matrix(row,col,nSamp).
%
% Examples:
%
%  Compute Noise Free and then multiple samples
%    sensorNF   = sensorComputeNoiseFree(sensor,oi);
%    voltImages = sensorComputeSamples(sensorNF,100);
%    imagesc(std(voltImages,0,3)); colorbar
%
%  To interact with a plane of the returned data you might use calls such
%  as
%    ii = 5
%    tmp = sensorSet(sensorNF,'volts',voltImages(:,:,ii));
%    e   = sensorGet(tmp,'electrons');
%
% Copyright ImagEval Consultants, LLC, 2011.

%% Define parameters
if notDefined('sensorNF'), error('Noise free sensor array required.'); end
if notDefined('nSamp'),  nSamp = 10; end
if notDefined('noiseFlag'), noiseFlag = 2; end  % Photon and electrical
if notDefined('showBar'), showBar = ieSessionGet('waitbar'); end

%%  Get noise free values
sz = sensorGet(sensorNF,'size');

%% Loop on number of samples to compute only the noise (no reuse)
sensorNF  = sensorSet(sensorNF,'noise flag',noiseFlag);  % 1 = photon, 2 = all
sensorNF  = sensorSet(sensorNF,'reuse noise',0); % Don't want to reuse

voltImages = zeros(sz(1),sz(2),nSamp);
str = sprintf('Computing %d samples',nSamp);
if showBar, h = waitbar(0,str); end

% Here, we use a simplied method to compute noise samples for human. The
% method is the same as the original one. But it doesn't do it in loop and
% checks are simplified. This could save a great amount of time in the
% whole simulation process
if sensorCheckHuman(sensorNF)
    error('For human sensor, use coneAbsorptions instead');
else    
    for kk=1:nSamp
        sensorN = sensorComputeNoise(sensorNF,[]);
        voltImages(:,:,kk) = sensorGet(sensorN,'volts');
        % v2 = voltImages(:,:,kk); v1 = sensorGet(sensorNF,'volts');
        % vcNewGraphWin; hist(v1(:)-v2(:),100)
        if ~mod(kk,10) && showBar, waitbar(kk/nSamp); end
    end
end
if showBar, close(h); end

end