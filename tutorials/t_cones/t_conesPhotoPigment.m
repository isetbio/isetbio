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
coneParams.indDiffParams.lambdaMaxShift = [-5 2 -0.5]';

% ISETBIo field of view.
%
% This determines the size of the scene and cone mosaic we work with, but
% not cone properties because as we will see below we are setting those
% explicitly in this script.
fovDegreesISETBio = 0.2;

% Compute the quantal sensitivity as probability of excitation given a
% photon of different wavelengths passing through the cone's entrance
% aperture.  Sensitivity is referred to srate of arrival of photons at the
% cornea, in that lens transmittance is folded in.
[~,T_quantalAbsorptionProb,T_quantalExcitationProb,adjIndDiffParams,cieConeParams,cieStaticParams] = ...
    ComputeCIEConeFundamentals(MakeItS(wls),coneParams.fieldSizeDegrees,...
    coneParams.ageYears,coneParams.pupilDiamMM,[],[],[], ...
    [],[],[],coneParams.indDiffParams);

% Plot
fundamentalsFig = figure; clf; hold on;
plot(wls,T_quantalExcitationProb(1,:),'r','LineWidth',6);
plot(wls,T_quantalExcitationProb(2,:),'g','LineWidth',6);
plot(wls,T_quantalExcitationProb(3,:),'b','LineWidth',6);
xlabel('Wavelength (nm)');
ylabel('Exication Probability');

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

% Numerical check. For some reason there is a small numerical difference
% between what we compute here and what was computed in the subroutine.
% Might have to do with order of operations and/or some difference
% between what Matlab does in a script versus a function. So here and below
% we don't check for exact equality but instead agreement to better than
% one part in 10^12, referenced to mean value.
if (max(abs(T_quantalAbsorptionProbFromPieces(:) - T_quantalAbsorptionProb(:)))/mean(T_quantalAbsorptionProb(:)) > 1e-12)
    error('Fail to compute quantal absorption probability same way twice');
end

% Add to plot and check numerically
figure(fundamentalsFig);
plot(wls,T_quantalExcitationProbFromPieces(1,:),'k','LineWidth',4);
plot(wls,T_quantalExcitationProbFromPieces(2,:),'k','LineWidth',4);
plot(wls,T_quantalExcitationProbFromPieces(3,:),'k','LineWidth',4);
if (max(abs(T_quantalExcitationProbFromPieces(:) - T_quantalExcitationProb(:)))/mean(T_quantalExcitationProb(:)) > 1e-12)
    error('Fail to compute quantal excitation probability same way twice');
end

%% Set up object for computing retinal image
%
% Turn off cos4 off axis scaling
oiBaseline = oiCreate('wvf human', coneParams.pupilDiamMM , [], wls);
opticsTemp = oiGet(oiBaseline,'optics');
opticsTemp = opticsSet(opticsTemp,'off axis method', 'Skip');
oiBaseline = oiSet(oiBaseline,'optics',opticsTemp);

% The Lens object in the oi object determines the lens transmittance used
% in the ISETBio computations. The code here shows how to set the lens
% density/transmittance to match the above.
%
% The spectral values in the unit density variable are multiplied by the
% scalar in the density variable to get the overall spectral density, and
% from thence the transmittance.  See "help Lens". This is
% overparameterized, but such parameterization is conventional so that the
% scalar density is separated from the absorbance (normalized to peak of
% 1).
% 
% Given the parameterization, the easiest way to adjust lens density is to
% set the unit density to the density we computed above and the density to
% 1. This makes a mockery of the convention that we separate absorbance
% from the scalar optical density.  One could normalize the density and use
% normalization factor to comptue the scalar optical density, and get the
% same answer.
lensObject = Lens('wave',wls,'unitDensity',-log10(lensTransmittance),'density',1);
oiBaseline = oiSet(oiBaseline,'lens',lensObject);

% Set up macular pigment object.  See comments where we set up lens object
% above.
macObject = Macular('wave',wls,'unitDensity',-log10(macTransmittance),'density',1);

% Set up photopigment object
%
% This is based on the axialDensity (parameter name 'opticalDensity'),
% absorbance (parameter name 'absorbance'), and quantalEfficiency
% (parameter name 'peakEfficiency'). Someday we will rewrite so that naming 
% conventions are consistent across different places we implement them, perhaps.
% 
% ISETBio wants absorbance in columns, while PTB wants it in rows, so we
% also need to transposes.
photopigmentObject = cPhotoPigment('wave', wls,...
    'opticalDensity',axialDensity,'absorbance',photopigmentAbsorbance', ...
    'peakEfficiency',cieStaticParams.quantalEfficiency );

%% Generate the mosaic object, with custom pigments
%
% This should now compute with the photoreceptor properties we defined
% above, which we will attempt to verify below.
cm = cMosaic(...
    'wave', wls, ...                % wavelength sampling
    'sizeDegs', [fovDegreesISETBio fovDegreesISETBio], ... % x,y size in degrees
    'eccentricityDegs', [0 0], ...  % eccentricity
    'tritanopicRadiusDegs', 0, ...  % Want some S cones and don't care about space in this example
    'macular', macObject, ...       % custom macular pigment object
    'pigment', photopigmentObject, ... % custom photopigment object
    'eccVaryingConeAperture', false, ... % for purposes here, keep cones the same size, etc.
    'eccVaryingConeBlur', false, ..., 
    'eccVaryingOuterSegmentLength', false, ...
    'eccVaryingMacularPigmentDensity', false, ...
    'eccVaryingMacularPigmentDensityDynamic', false, ...
    'anchorAllEccVaryingParamsToTheirFovealValues', true ...
    );
LConeIndices = find(cm.coneTypes == 1);
MConeIndices = find(cm.coneTypes == 2);
SConeIndices = find(cm.coneTypes == 3);

% Obtain ISETBio's view of cone aperture size
coneDiameterListUM = cm.coneApertureDiametersMicrons;
if any(diff(coneDiameterListUM) ~= 0)
    error('cm cone diameters do not all have the same size');
end
coneApertureDiamUM = coneDiameterListUM(1);
coneAreaUM2 = pi*(coneApertureDiamUM/2)^2;

% Get the ISTBio integration time
integrationTimeSec = cm.integrationTime;

%% Generate the baseline stimuli
%
% A set of monochromatic images at each sample wavelength.  The photon
% level is scene radiance in photons/sr-m2-nm-sec.
pixelSize = 64;
photonsPerSrM2NMSec = 1e25;
for ww = 1:length(wls)
    % Set up a dummy scene. Spatially uniform with a black body
    % spectrum of 5000 degK.  We'll replace the scene contents just below.
    scene{ww} = sceneCreate('uniform bb',pixelSize,5000,wls);

    % Use small field of view to minimize effects of eccentricity, and also
    % so we don't need too many pixels (for efficiency in this demo).
    scene{ww} = sceneSet(scene{ww},'fov',fovDegreesISETBio);

    % Get photons and rewrite to be monochromatic constant power in
    % photons/sec-nm.
    photons = sceneGet(scene{ww},'photons');
    photons = zeros(size(photons));
    photons(:,:,ww) = photonsPerSrM2NMSec*ones(pixelSize,pixelSize);
    scene{ww} = sceneSet(scene{ww},'photons',photons);
end

%% Compute the retinal image and cone excitations for each scene
%
% Use this to get the cone fundamental that ISETBio is effectively
% using.
T_quantalExcitationISETBio = zeros(size(T_quantalExcitationProb));
for ww = 1:length(wls)
    % Compute retinal image
    oiComputed{ww} = oiCompute(scene{ww}, oiBaseline);

    % Compute noise free cone excitations
    coneExcitations{ww} = cm.compute(oiComputed{ww},'nTrials', 1);

    % Find L, M and S cone excitations at this wavelength.  This is
    % accomplished by extracting the mean response of each cone type from
    % the excitations computed just above and diviting by the input power
    % in the scene.
    %
    % We expect the answer to be proportional to the fundamentals we
    % computed above, because we have not (yet) accounted for the geometry
    % between radiance in the scene and retinal irradiance, nor the cone
    % integration time.
    excitationsConeISETBioMono(1,ww) = mean(coneExcitations{ww}(LConeIndices));
    excitationsConeISETBioMono(2,ww) = mean(coneExcitations{ww}(MConeIndices));
    excitationsConeISETBioMono(3,ww) = mean(coneExcitations{ww}(SConeIndices));
    T_quantalExcitationISETBio(1,ww) = excitationsConeISETBioMono(1,ww)/photonsPerSrM2NMSec;
    T_quantalExcitationISETBio(2,ww) = excitationsConeISETBioMono(2,ww)/photonsPerSrM2NMSec;
    T_quantalExcitationISETBio(3,ww) = excitationsConeISETBioMono(3,ww)/photonsPerSrM2NMSec;
end
excitationsConeISETBio = sum(excitationsConeISETBioMono,2);

% Scale and plot what we get this way
scaleFactor = T_quantalExcitationISETBio(:)\T_quantalExcitationProb(:);
figure(fundamentalsFig);
plot(wls,scaleFactor*T_quantalExcitationISETBio(1,:),'y-','LineWidth',3);
plot(wls,scaleFactor*T_quantalExcitationISETBio(2,:),'y-','LineWidth',3);
plot(wls,scaleFactor*T_quantalExcitationISETBio(3,:),'y-','LineWidth',3);

% This also matches quantum efficiency in the cm object, with lens
% added back in.
T_quantalExcitationISETBioCM = cm.qe';
plot(wls,T_quantalExcitationISETBioCM(1,:) .* lensTransmittance,'b:','LineWidth',1.5);
plot(wls,T_quantalExcitationISETBioCM(2,:) .* lensTransmittance,'b:','LineWidth',1.5);
plot(wls,T_quantalExcitationISETBioCM(3,:) .* lensTransmittance,'b:','LineWidth',1.5);

%% Get parameters needed for retinal irradiance from scene radiance in PTB
%
% Start by getting parameters we need from the optics.
opticsBaseline = oiGet(oiBaseline,'optics');
if (opticsGet(opticsBaseline,'pupil diameter','mm') ~= coneParams.pupilDiamMM)
    error('Baseline optics does not have expected pupil diameter');
end
focalLengthMM = opticsGet(opticsBaseline,'focal length','mm');
pupilAreaMM2 = pi*(coneParams.pupilDiamMM/2)^2;

% Also get rradiance correction factor between PTB and ISETBio
%
% See function for oiCalculateIrradiance for where this happens in ISETBio
% and some info on why. This correctrion applied in that routine remains
% mysterious to some of us, but for our purposes here we just match the calculation.
%
% This is what ISETBio does:
%     irradiance = pi / (1 + 4 * fN ^ 2 * (1 + abs(m)) ^ 2) * radiance;
% This is what PTB does:
%     irradiance = pupilArea/(eyeLength^2) * radiance;
%                = (pi/4)*(pupilDiam/eyeLength)^2 * radiance
%                = pi/(4 * fN ^ 2)
% because fN = eyeLength/pupilDiam.
fN = focalLengthMM/coneParams.pupilDiamMM;
m = opticsGet(opticsBaseline,'magnification',sceneGet(scene{1},'distance'));
irradianceCorrectionFactorPTBToISETBio = (4*fN^2)/(1 + 4 * fN^2 * (1 + abs(m)) ^ 2);

%% Convert radiance of source to retinal irradiance and then to cone excitations.
%
% excitations, using PTB method.  Do this at each wavelength, so that we
% can then compare to the same calculation done in ISETBio above.
excitationsConePTBMono = zeros(size(excitationsConeISETBioMono));
for ww = 1:length(wls)

    % Set up monochromatic radiance and convert to irradiance in ISETBio
    radiancePhotonsPerSecM2Sr = zeros(size(wls));
    radiancePhotonsPerSecM2Sr(ww) = photonsPerSrM2NMSec;
    irradiancePhotonsPerSecUM2 = RadianceToRetIrradiance(radiancePhotonsPerSecM2Sr,wls, ...
        pupilAreaMM2,focalLengthMM);

    % Correct irradiance to match what ISETBio does.
    irradiancePhotonsPerSecUM2Corrected = irradiancePhotonsPerSecUM2*irradianceCorrectionFactorPTBToISETBio;
    irradiancePhotonsPerSecM2CorrectedChk(ww) = sum(irradiancePhotonsPerSecUM2Corrected .* lensTransmittance' * 1e12);

    % Get cone exciations from retinal irradiance.  The approximation to
    % the integral here more or less matches the way it is done in ISETBio. 
    for cc = 1:3
        excitationsConeSec(cc) = T_quantalExcitationProb(cc,:)*irradiancePhotonsPerSecUM2Corrected*coneAreaUM2*(wls(2)-wls(1));
        excitationsConePTBMono(cc,ww) = excitationsConeSec(cc)*integrationTimeSec;
    end
end

%% Plot and compare with ISETBio
excitationsFig = figure; clf; hold on;
plot(wls,excitationsConePTBMono(1,:),'r','LineWidth',6);
plot(wls,excitationsConePTBMono(2,:),'g','LineWidth',6);
plot(wls,excitationsConePTBMono(3,:),'b','LineWidth',6);
plot(wls,excitationsConeISETBioMono(1,:),'k','LineWidth',2);
plot(wls,excitationsConeISETBioMono(2,:),'k','LineWidth',2);
plot(wls,excitationsConeISETBioMono(3,:),'k','LineWidth',2);
xlabel('Wavelength (nm)');
ylabel('Excitations/Cone-Integration Time');
if (max(abs(excitationsConePTBMono(:) - excitationsConeISETBioMono(:))) ...
        /mean(excitationsConePTBMono(:)) > 1e-6)
    error('Fail to compute excitations same way in ISETBio and PTB, after correction');
end

% Compare PTB with ISETBio irradiance calculation
%
% Be sure to convert from UM2 to M2 and take lens transmittance into
% account.  The version of the PTB calculation here corrects for the known
% difference in how ISETBio does this calculation.
irradianceFig = figure; clf; hold on;
plot(wls,irradiancePhotonsPerSecM2CorrectedChk ,'k','LineWidth',4);
irradiancePhotonsPerSecM2ISETBio = zeros(size(irradiancePhotonsPerSecUM2Corrected));
for ww = 1:length(wls)
    temp = oiGet(oiComputed{ww},'photons');
    irradiancePhotonsPerSecM2ISETBio(ww) = mean(mean(temp(:,:,ww)));
    plot(wls(ww),irradiancePhotonsPerSecM2ISETBio(ww),'ro','MarkerFaceColor','r','MarkerSize',10);
end
xlabel('Wavelength (nm)'); ylabel('Retinal Irradiance (photons/sec-m2)');
if (max(abs(irradiancePhotonsPerSecM2CorrectedChk(:) - irradiancePhotonsPerSecM2ISETBio(:))) ...
        /mean(irradiancePhotonsPerSecM2ISETBio(:)) > 1e-6)
    error('Fail to compute retinal irradiance same way in ISETBio and PTB, after correction');
end

    
