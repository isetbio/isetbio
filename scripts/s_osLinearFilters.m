%% Calculate he lms linear filters in the BioPhys model (Rieke)
%
% The linear filters are calculated by introducing a small perturbation
% into the osBioPhys model.  The linear filters differ as a function of the
% mean luminance level (background absorption rate).
%
% Notice in the resulting
%
% Computation - We create a uniform scene and set its luminance level.
% Then we create a cone mosaic of one type of 
%

%%
ieInit

%%
scene = sceneCreate('uniform ee');
scene = sceneSet(scene,'fov',1);
oi = oiCreate('human');

%%
cMosaic = coneMosaic('os',osLinear,'spatialDensity',[0 1 0 0]);
cMosaic.setSizeToFOV(0.7);
cMosaic.integrationTime = 0.001;

% Initialize luminance levels and mean isomerization count and filters
lum = logspace(log10(0.5),4,12);
meanI = zeros(size(lum));
f = zeros(4002,length(lum));

%%  Loop to calculate filters

for ii = 1:length(lum)
    
    scene = sceneSet(scene,'mean luminance',lum(ii));
    oi = oiCompute(oi,scene);
    % oiGet(oi,'mean illuminance')
    
    cMosaic.compute(oi);
    % mean(cMosaic.absorptions(:))   
    tmp = coneMeanIsomerizations(cMosaic);
    meanI(ii) = tmp(1);    % Mean number of isomerizations
    
    % Compute the linear filters and plot them
    cMosaic.os.linearFilters(cMosaic);
    f(:,ii) = cMosaic.os.lmsConeFilter(:,1);
    
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

loglog(meanI,1./mx,'k-o');            
xlabel('Mean isomerization rate');
ylabel('Effective threshold');
grid on;
title('Threshold vs. Background Rate')

%%


