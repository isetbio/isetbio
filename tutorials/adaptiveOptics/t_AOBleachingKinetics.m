%% Bleaching kinetics in ISETBio
%
% Description:
%    This reproduces Figure S1 (panel D) of Cooper et al. (2020), showing how the
%    fraction of L/M cone pigment bleaches evolves over time with a
%    stimulus train of 1 sec flashes.
%
% See also: t_AODisplay, t_radianceToCornealIrradiance

% History:
%   01/09/21  dhb  Wrote it.

%% Initialize
close all; clear; ieInit;

%% Stimulus spectral properties
%
% Narrowband stimulus modeled as Gaussians centered on specified
% wavelength, with specified FWHM
whichStimCondition = 'stimulus';
switch (whichStimCondition)
    case 'stimulus'
        stimulusWl = 545;                         % Wavelength of narrowband background field. This can be the imaging light
        stimulusCornealPowerNw = 450;             % Stimulus power passing through pupil in nanowatts
        stimulusCornealPowerUw = stimulusCornealPowerNw*(1e-3);
        
    case 'background'
        stimulusWl = 790;                         % Wavelength of narrowband background field. This can be the imaging light
        stimulusCornealPowerUw = 111;             % Stimulus power passing through pupil in microwatts
    
    otherwise
        error('Unknown stimulus condition specified');
end
stimulusFWHM = 5;                         % Full-width at half max of background (in nm).
 
%% Wavelength support
deltaWl = 1;
wls = (500:deltaWl:900)';
nWls = length(wls);

%% Stimulus spatial parameters
backgroundLinearSizeDegs = 1;              % Linear side of square field.

%% Make relative spectral power distributions.
%
% Approximated by a Gaussian with specified center wavlength and FWHM,
% and with total power given by the corneal power specified above.
% The call to trapz takes wavelength spacing into account when normalizing the power.
relSpd = normpdf(wls,stimulusWl,FWHMToStd(stimulusFWHM));
unitSpd = relSpd/trapz(wls,relSpd);

%% Pupil size
pupilDiameterMm = 7;                        % Pupil diameter.
pupilAreaMm = pi*(pupilDiameterMm/2)^2;     % For convenience below, compute pupil area.

%% Convert stimulus power to trolands
%
% Get equivalent spectral radiance as a function of the wavelength support.
%
% The routine here finds the radiance on an external conventional display that
% produces the same retinal illuminance as the corneal power specified
% above.  This is purely geometric calculation; attenuation of light by
% occular media is not taken into account at this stage.  Note that this
% conversion routine expects power per wavelength band, not power per nm,
% as its input, but returns power per nm as its output.  A little
% confusing, but that is why the spd being passed in is multiplied by the
% wavelength spacing.
stimulusSpdCornealPowerUw = stimulusCornealPowerUw*unitSpd;
stimulusSpdRadiance = AOMonochromaticCornealPowerToRadiance(wls,wls,stimulusSpdCornealPowerUw*deltaWl,pupilDiameterMm,backgroundLinearSizeDegs^2);

% Make sure our computed radiance yields the desired corneal
% irradiance when we go in the other direction.  The magic
% numbers (1e6) and (1e-3) in the call just below do unit conversions
% from units we're using here to those expected by RadianceAndDegrees2ToCornIrradiance
stimulusSpdCornealIrradianceUWMm2Check = RadianceAndDegrees2ToCornIrradiance(stimulusSpdRadiance,backgroundLinearSizeDegs^2)*(1e6)*((1e-3)^2);
stimulusCornealPowerUWCheck = trapz(wls,stimulusSpdCornealIrradianceUWMm2Check)*pupilAreaMm;
if (abs(stimulusCornealPowerUWCheck-stimulusCornealPowerUw)/stimulusCornealPowerUw > 1e-4)
    error('Do not get right cornal power back from computed radiance');
end

%% Convert from radiance to luminance.
%
% The magic number 683 makes the unis of luminance cd/m2, given the
% Watts/[sr-m2-nm] units we're using for radiance.
load T_xyz1931
T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,wls);
T_y = T_xyz(2,:);
stimulusLuminance = trapz(wls,T_y'.*stimulusSpdRadiance);
fprintf('Stimulus luminance is %0.1f cd/m2\n',stimulusLuminance);

%% Convert luminance to trolands for bleaching calcs
stimulusTrolands = stimulusLuminance*pupilAreaMm;
fprintf('Stimulus is %0.1f trolands\n',stimulusTrolands);

% Compute a temporal stimulus timecourse
timePerTrialSec = 18;
stimStartTimeSec = 4;
timePerStimSec = 1;
nTrials = 13;
totalTimeSec = timePerTrialSec*nTrials;
totalTimeMsec = 1000*totalTimeSec;
timeMsec = (1:totalTimeMsec) - 1;
timeSec = timeMsec/1000;

switch (whichStimCondition)
    case 'stimulus'
        trolands = zeros(size(timeMsec));
        for ii = 1:nTrials
            startTime = (ii-1)*timePerTrialSec*1000 + stimStartTimeSec*1000;
            finishTime = startTime + timePerStimSec*1000;
            index = find(timeMsec > startTime & timeMsec < finishTime);
            trolands(index) = stimulusTrolands;
        end
        
    case 'background'
        trolands = stimulusTrolands*ones(size(timeMsec));
                
    otherwise
        error('Unknown stimulus condition specified');
end


%% Compute bleaching over time
initialFractionBleached = 0;
fractionBleached = ComputePhotopigmentBleaching(trolands,'cones','trolands','Boynton',initialFractionBleached,'msec');
fractionUnbleached = 1 - fractionBleached;

%% Plot of bleaching
figure; clf;
subplot(2,1,1); hold on;
plot(timeSec,trolands,'g','LineWidth',4);
xlim([0 totalTimeSec]);
xlabel('Time (sec)');
ylabel('Trolands');
subplot(2,1,2); hold on;
plot(timeSec,fractionUnbleached,'g','LineWidth',4);
xlim([0 totalTimeSec]);
ylim([0 1]);
xlabel('Time (sec)');
ylabel('Fraction L/M cone pigment unbleached');