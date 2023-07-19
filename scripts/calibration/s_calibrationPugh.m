%% Calibration testing
%
% For unit testing
%
%  Ed Pugh says (Pugh_Summary_for_FFB.pdf) that based on Delori and Webb
%  2007 he knows that
%
% At the eye:
%   10^15 photons cm^-2 sec^-1 at 580 nm corresponds to a retinal
%   illuminance of 590,000 photopic Trolands.
%
% He then says that a photopic Troland is related to candelas/m2 by the
% idea that a 1 cd m^-2 patch seen a 1 mm^2 pupil area.  For a pupil that
% is 2 mm diameter, the area is 3.1mm^2.  So converting back to the scene,
%
%   590,000 Trolands at the retina at 580nm is supposed to be 190,000
%   cd/m2 at the scene
%
% Let's check these numbers by
%
%  1. Creating a 580 nm uniform scene 
%  2. Imaging it through a lens with a 1mm pupil, focal length of 17mm
%  3. Calculating the retinal illuminance in the optical image window
%
% We should find that
%     190,000 cd/m2 at 580 nm in the scene, produces
%     10^19 photons m^-2 s^-1 at the retina.
%
%%
ieInit;
 
meanluminance   = 190000;  % Cd/m2
monochromeWave  = 580;     % 
 
%% set uniform monochromatic scene
 
scene = sceneCreate('uniformmonochromatic');
scene = initDefaultSpectrum(scene,'custom', monochromeWave);
scene = sceneSet(scene,'mean luminance', meanluminance);   % Cd/m2
sceneWindow(scene);
 
%% create an optical image of human eye
oi = oiCreate('human');
optics = opticsCreate('human', 1e-3);  % 1 mm
oi = oiSet(oi,'optics',optics);
 
oi = oiCompute(oi,scene);
vcAddAndSelectObject(oi);
oiWindow(oi);
 
%% calc irradiance at 580 nm
centerPoints = round(size(oi.depthMap) ./ 2 ) ;
irradiance = vcGetROIData(oi,centerPoints,'photons');
irradiance = mean(irradiance);  % quanta / sec / m^2 /nm

%% Print
fprintf('Should be close to 10^19\nIrradiance: %.3e (q/s/m^2/nm)  at %.0f nm\n\n',irradiance,monochromeWave);
fprintf('Should be close to 10^19\nIrradiance: 10^%1.1f (q/s/m^2/nm)  at %.0f nm\n\n',...
    log10(irradiance),monochromeWave);
 
%%

