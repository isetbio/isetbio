% Illustrate photoPigment object
%
% Description:
%    The photoPigment object represents the spectral responsivity of
%    the cone photopigments.  These values are needed to calculate the
%    capture of light by cones, along with the oi spectral irradiance
%    and the other inert pigments (lens, macular). 
%
%    This tutorial illustrates how it works.
%
% See also
%    

%% Initialize
ieInit; clear; close all;

%% Set wavelength support
wls = (400:5:700)';

%% Start by computing outside of ISETBio
%
% With this call, we get cone fundamentals in quantal units, so
% that integrating against spectral retinal irradiance and multiplying
% by cone aperture acceptance area gives excitations/cone-sec. Note
% that these take lens density into account.
%
% The fundamentals are according to the CIE standard but also include the
% Asano et al. individual differences model.
%
% What we'll do below is show how to use the various pieces available from
% this call below, plus the cone aperture size, inside of ISETBio.
coneParams = DefaultConeParams('cie_asano');

% Set age, pupil size, and field size paramters.  These are part of the CIE
% model, and we apply the Asano et al. adjustments to these.
coneParams.fieldSizeDegrees = 2;
coneParams.ageYears = 32;
coneParams.pupilDiamMM = 3;

% Set individual difference parameters. These are the parameters we
% implement from the Asano et al. model. Density changes apply to lens,
% macular pigment, and L, M, S photopigment. These are expressed as percent
% changes from the nominal values. Shifts in photoreceptor wavelength of
% peak absorption are in nm.
%
% These are set to non-zero here just to make sure tests of agreement below
% aren't specific to the more generic case of no shifts.
coneParams.indDiffParams.dlens = 1;
coneParams.indDiffParams.dmac = -2;
coneParams.indDiffParams.dphotopigment = [5 -3 6]';
coneParams.indDiffParams.lambdaMaxShift = [-0.5 2 -10]';

% Compute the quantal sensitivity as probability of excitation given a
% photon of different wavelengths passing through the cone's entrance
% aperture.  Sensitivity is referred to srate of arrival of photons at the
% cornea, in that lens transmittance is folded in.
[~,T_quantalAbsorptionProb,T_quantalExcitationProb,adjIndDiffParams,cieConeParams,cieStaticParams] = ...
    ComputeCIEConeFundamentals(MakeItS(wls),coneParams.fieldSizeDegrees,...
    coneParams.ageYears,coneParams.pupilDiamMM,[],[],[], ...
    [],[],[],coneParams.indDiffParams);

%% Obtain ISETBio's view of cone aperture size at the fovea
%
% See coneDensityReadData for explanation of retinal coordinate system.
eccentricity = 0.0;
angle = 0;
[~, coneApertureDiamM, ~] = coneSizeReadData('eccentricity', eccentricity, 'angle', angle);
coneApertureDiamUm = coneApertureDiamM*1e6;
coneAreaUm2 = pi*(coneApertureDiamUm/2)^2;

%% Get fundamentals taking cone aperture into account
%
% These express the cone excitation rate when retinal irradiance is in units
% of photons/um^2-sec, but expressed at the cornea (i.e. not including lens
% transmittance).
%
% Below we'll build these up from parts to illustrate the calculations in
% more detail.
T_quantalExcitationIrradiance = T_quantalExcitationProb*coneAreaUm2;

% Plot
figure; clf; hold on;
plot(wls,T_quantalExcitationIrradiance(1,:),'r','LineWidth',4);
plot(wls,T_quantalExcitationIrradiance(2,:),'g','LineWidth',4);
plot(wls,T_quantalExcitationIrradiance(3,:),'b','LineWidth',4);
xlabel('Wavelength (nm)');
ylabel('Fundamental');

%% Lens transmittance
%
% Get lens transmittance.  This call computes lens density according to the
% CIE standard, adjusting for age and pupil diameter, and then returns
% transmittance as 10.^(-lensDensity).
% 
% The calculation builds up the density from two tabulated parts, one of
% which is adjusted for age through multiplication with an age dependent
% factor, and the other which is age independent.  The result is then
% adjusted multiplicativey by a pupil diameter dependent factor.
lensTransmittance = LensTransmittance(wls,'Human','CIE', ...
    coneParams.ageYears,coneParams.pupilDiamMM);

% Adjust lens density according to the individual differenes
% Asano model. This specifies a percentage change from the default lens
% density for a given age and pupil size.  The calculation works
% by going from transmittance back to density, applying the specified
% factor, and then converting back to
% transmittance.
unadjustedLensDensity = -log10(lensTransmittance);
adjustedLensDensity = unadjustedLensDensity * (1 + coneParams.indDiffParams.dlens/100);
lensTransmittance = 10.^-adjustedLensDensity;

% Check that our computation here matches what was done inside of
% ComputeCIEConeFundamentals.
if (any(lensTransmittance ~= adjIndDiffParams.lens))
    error('Fail to compute lens transmittance the same way twice');
end

%% Macular pigment transmittance
%
% The exercise parallels that for lens density above and is not broken out
% or commented as fully.
macTransmittance = MacularTransmittance(wls,'Human','CIE',coneParams.fieldSizeDegrees);
unadjustedMacDensity = -log10(macTransmittance);
adjustedMacDensity = unadjustedMacDensity * (1 + coneParams.indDiffParams.dmac/100);
macTransmittance = 10.^-adjustedMacDensity;
if (any(macTransmittance ~= adjIndDiffParams.mac))
    error('Fail to compute macular transmittance the same way twice');
end

%% Photopigment absorptance
%
% Start with the absorbance, which is density normalized to a peak of 1.
% This is a tabulated function for each cone photopigment type, with the
% tables holding log10 absorbance.
temp = load('T_log10coneabsorbance_ss');
photopigmentAbsorbance = 10.^SplineCmf(temp.S_log10coneabsorbance_ss,temp.T_log10coneabsorbance_ss,wls,2);
clear temp

% The individual differences model allows for shifting the absorbance along
% the wavelength axis. Do that here.
photopigmentAbsorbance = ShiftPhotopigmentAbsorbance(wls,photopigmentAbsorbance, ...
    coneParams.indDiffParams.lambdaMaxShift,coneParams.indDiffParams.shiftType);
if (any(adjIndDiffParams.absorbance(:) ~= photopigmentAbsorbance(:)))
    error('Fail to compute photopigment absorbance same way twice');
end

% Get axial density.  This scales the normalized photopigment absorbance we
% computed just above.
axialDensity = cieConeParams.axialDensity .* (1 + coneParams.indDiffParams.dphotopigment/100);
if (any(adjIndDiffParams.dphotopigment(:) ~= axialDensity(:)))
    error('Mysterious axial density difference');
end

% Compute absorbtance from absorbance.  The absorbtance is the probabilyt that a photon
% entering the cone will be absorbed, so 1- tranmittance = 1 - 10^(-absorbance*axialDensity).
% This calculation is done by routine AbsorbanceToAbsorbtance;
photopigmentAbsorptance = AbsorbanceToAbsorptance(photopigmentAbsorbance,wls,axialDensity);
if (any(adjIndDiffParams.absorptance(:) ~= photopigmentAbsorptance(:)))
    error('Fail to compute photopigment absorptance same way twice');
end

%% Put it all together
%
% We multiply the lens transmittance, the macular transmittance, and the
% photopigment absorbtance to get the probability that a photon starting at
% the cornear and traced into a cone aperture will be absorbed.
%
% Not all absorptions lead to an excitation (aka isomerization of a
% photopigment molecule.) The probability an absorption does lead to an
% excitation is the quantum efficiency.  We incorporate that here too.
for cc = 1:3
    T_quantalAbsorptionProbFromPieces(cc,:) = lensTransmittance .* macTransmittance .* photopigmentAbsorptance(cc,:);
    T_quantalExcitationProbFromPieces(cc,:) = T_quantalAbsorptionProbFromPieces(cc,:) * cieStaticParams.quantalEfficiency(cc);
end
if (max(abs(T_quantalAbsorptionProbFromPieces(:) - T_quantalAbsorptionProb(:)))/mean(T_quantalAbsorptionProb(:)) > 1e-12)
    error('Fail to compute quantal absorption probability same way twice');
end

% Irradiance version
T_quantalExcitationIrradianceFromPieces = T_quantalExcitationProbFromPieces*coneAreaUm2;

% Add to plot and check numerically
plot(wls,T_quantalExcitationIrradianceFromPieces(1,:),'k','LineWidth',2);
plot(wls,T_quantalExcitationIrradianceFromPieces(2,:),'k','LineWidth',2);
plot(wls,T_quantalExcitationIrradianceFromPieces(3,:),'k','LineWidth',2);
if (max(abs(T_quantalExcitationProbFromPieces(:) - T_quantalExcitationProb(:)))/mean(T_quantalExcitationProb(:)) > 1e-12)
    error('Fail to compute quantal excitation probability same way twice');
end

%% Generate the ring rays stimulus
radialFrequency = 8;
pixelSize = 256;
scene = sceneCreate('rings rays',radialFrequency,pixelSize,wls);
scene = sceneSet(scene, 'fov', 0.5);

%% Compute the retinal image
%
% First create the optical image object.
pupilDiamterMm = 3;
oi = oiCreate('wvf human', pupilDiamterMm , [], wls);

% Compute the optical image
oi = oiCompute(scene, oi);

% The lens density in the oi object affect the cone fundamentals.  This
% code shows how to adjust lens density.

% The values in the unit density variable are multiplied by the scalar in
% the density variable to get the overall spectral density, and from thence the
% transmittance.  See "help Lens".
lens0 = oiGet(oi,'lens');
lensUnitDensity1 = lens0.unitDensity;
lensPeakDensity1 = lens0.density;
lens1 = Lens('wave',wls,'unitDensity',lensUnitDensity1,'density',lensPeakDensity1);
oi = oiSet(oi,'lens',lens1);

%% Generate the mosaic
cm = cMosaic(...
    'wave', wls, ...                % wavelength sampling
    'sizeDegs', [0.5 0.5], ...      % x,y size in degrees
    'eccentricityDegs', [0 0], ...  % eccentricity
    'eccVaryingConeBlur', false ... % for purposes here, keep cones the same size
    );

%% Create photoPigment object
%
% Plot spacing at 5 nm looks a little nicer than the default spacing.
pp = photoPigment('wave', wls);

%% Absorbance
% The cone absorbance function used by default is obtained via data routine
% coneAbsorbanceReadData.
%
% Absorbance is sometimes called optical density.
%
% The peak absorbance is 1 by convention in how the values are tabulated.
ieNewGraphWin;
plot(pp.wave, pp.absorbance)
grid on;
xlabel('Wavelength (nm)');
ylabel('Relative sensitivity');



%%