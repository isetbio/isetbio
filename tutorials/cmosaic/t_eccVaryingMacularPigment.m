% t_eccVaryingMacularPigment
%
% This tutorial shows how to use our machinery
% to compute the macular pigment spectral transmittance
% for any eccentricity.

% History
%  11/05/24  dhb  Wrote it.

%% Initialize
clear; close all;

%% Decide wavelength spacing to be used (nm)
wls = (380:5:720)';

%% Set up macular pigment object.
%
% This is set up for standard foveal values
mac = Macular('wave',wls);

%% Plot foveal macular density and transmittance
%
% The actual foveal spectral density is given by the
% unitdensity times the density.  The unitDensity property is
% the spectra density normalized to a peak of 1 (before
% sampling at the discrete wavelengths specified) and
% the density property is the peak density, here defaulted
% to the foveal value.
figure; clf;
subplot(1,2,1); hold on;
plot(mac.wave,mac.unitDensity*mac.density,'r','LineWidth',3);
xlabel('Wavelength (nm)');
ylabel('Macular Pigment Density');
subplot(1,2,2); hold on;
plot(mac.wave,mac.transmittance,'r','LineWidth',3);
xlabel('Wavelength (nm)');
ylabel('Macular Pigment Transmittance');

%% Plot how peak density varies with eccentricity
% 
% The eccDensity method returns density as a function
% of eccentricity in degrees, assumed circularly symmetric.
% See "help Macular.eccDensity" for information on where
% the data came from and how it was smoothed.
%
% This is hard coded to a peak density of 0.35, but you
% can scale it if you like.
%
% Note that we do not try to force consistency between this
% estimate of how macular pigment density varies with eccentricity
% and the CIE cone fundamental standard as a function of eccentricity,
% Because 0.35 peak density matches the the CIE 2 deg fundamentals, it
% may be a little low at the very center.
eccentricitiesDegs = (0:0.1:10)';
figure; clf; hold on;
plot(eccentricitiesDegs,mac.eccDensity(eccentricitiesDegs),'r','LineWidth',3);
xlabel('Eccentricity (degs)');
ylabel('Macular pigment peak density');

%% Adust the eccentricity varying function for the peak density we are using.
peakDensity = mac.density;
eccPeakDensity = 0.35;
eccVaryingMacularDensity = (peakDensity/eccPeakDensity)*mac.eccDensity(eccentricitiesDegs);

%% Compute macular pigment transmittance at each eccentricity
macTransmittance = zeros(length(wls),length(eccentricitiesDegs));
for ee = 1:length(eccentricitiesDegs)
    % Direct calculation
    macTransmittance(:,ee) = 10.^-(eccVaryingMacularDensity(ee) * mac.unitDensity);

    % Check using macular object machinery
    mac.density = eccVaryingMacularDensity(ee);
    checkTransmittance = mac.transmittance;
    if (max(abs(macTransmittance(:,ee)-checkTransmittance)) > 1e-10)
        error('Do not understand macular pigment transmittance calculation');
    end
end

%% Plot macular pigment transmittance a a few chosen eccentricities
%
% Plot code is a little fragile if you start adding or changing
% eccentricities.
plotEccentricities = [0 2 4 6 8];
plotColors = ['r', 'c', 'k', 'g', 'b'];
theLegend = {};
figure; clf; hold on;
for ee = 1:length(plotEccentricities)
    index = find(eccentricitiesDegs == plotEccentricities(ee));
    plot(wls,macTransmittance(:,index),plotColors(ee),'LineWidth',3);
    theLegend{ee} = [num2str(plotEccentricities(ee)) ' deg'];
end
xlabel('Wavelength (nm)');
ylabel('Macular transmittance');
legend(theLegend);