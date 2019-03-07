%% The cone mosaic and the lms linear filters in the BioPhys model (Rieke)
%
% The linear outer segment model relies on impulse response functions that
% are calculated by introducing a small perturbation into the osBioPhys
% model.  The linear filters differ as a function of the mean luminance
% level (background absorption rate).
%
% This script illustrates the filters, their dependence on background rate,
% and how that dependence translates into a 'gain' (threshold sensitivity).
%
% This script carries out a very simple but relatively complete computation
% that uses the filters.
% 
%   *  We create a uniform scene and set its luminance level.
%   *  We create a cone mosaic of one type of cone (M)
%   *  For many different luminances, we calculate the linear filters
%   *  We plot the threshold as a function of luminance
%

%%
ieInit

%%
scene = sceneCreate('uniform ee');
scene = sceneSet(scene,'fov',1);
oi = oiCreate('human');

%%
cMosaic = coneMosaic('os',osLinear,'spatialDensity',[0 0 1 0]);
cMosaic.setSizeToFOV(0.7);
cMosaic.integrationTime = 0.001;

% Initialize luminance levels and mean isomerization count and filters
lum = logspace(log10(0.5),4,12);
meanI = zeros(size(lum));
tSamples = size(cMosaic.os.linearFilters(cMosaic),1);

f = zeros(tSamples,length(lum));

%%  Loop to calculate filters

for ii = 1:length(lum)
    
    scene = sceneSet(scene,'mean luminance',lum(ii));
    oi = oiCompute(oi,scene);
    % oiGet(oi,'mean illuminance')
    
    cMosaic.compute(oi);
    % mean(cMosaic.absorptions(:))   
    tmp = coneMeanIsomerizations(cMosaic);
    meanI(ii) = tmp(2);    % Mean number of isomerizations
    
    % Compute the linear filters and plot them
    cMosaic.os.linearFilters(cMosaic);
    f(:,ii) = cMosaic.os.lmsConeFilter(:,2);
    
end
%% Do the plot

vcNewGraphWin([],'wide');
hold on;
t = cMosaic.os.timeAxis;

subplot(1,2,1)
plot(t,f);
grid on
xlabel('Time (sec)'); ylabel('Current (pA)');

l = cell(1,length(lum));
for ii=1:length(meanI), l{ii} = num2str(meanI(ii),'%.1e'); end
legend(l)

title('Impulse Response vs. Background Rate')

%% Calculate the peak of the impulse response

subplot(1,2,2)
mx = max(f);  

% The peak response is basically the impulse response gain. To be at
% threshold, the photons times the gain equals some criterion level
%
%   ThreshPhotons * Gain = Criterion signal
%
% So 
%
%   ThresholdPhotons = Criterion / Gain
%

loglog(lum,1./mx,'k-o');            
xlabel('Mean luminance (cd/m^2)');
ylabel('Effective threshold');
grid on;
title('Threshold vs. Background Rate')

%%


