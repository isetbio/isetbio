%% s_calibrationRough
%
% Calculate the cone absorptions rates per second to a range of blackbody
% radiators.
%
% Where does the rule of thumb below come from? Someone seems to have
% told this to Wandell at some point in the past.
%
% A 100 cd/m2 5000 deg blackbody radiator causes about 10,000
% absorptions in a second in the L cones according to someone in BW's
% past.  
% 
% For the very small foveal cones in this simulation the number is
% about 7900, and if you move about 0.5 deg to the periphery the
% number becomes about 20,000.  So BW views the 10K as a good
% approximation.
%
% A second rule of thumb BW once heards is that the L,M,S abosrptions
% are 95,66,7 absorptions/sec per candela/m2 for the L, M and S cones.
%
% This is within the range as well, for the cone aperture sizes near
% the fovea and just outside.
%
% Maybe DHB knows where BW heard these numbers.  They are obviously
% just rough values, as represented in the name of this script.
%
% BW Vistasoft Team, 2014

%%
ieInit;

%%
lum       = 100;                           % cd/m2
noiseFlag = 0;

oi = oiCreate('human');

% The aperture size, and thus number of absorptions, varies with
% eccentricity.
cm = coneMosaicRect('eccentricity units','deg','center',[0,0]);
cm.setSizeToFOV(2);      % Small FOV
cm.integrationTime = 1;  % One second

%% Set variables
cTemp = 5000;   % Color temperature
scene = sceneCreate('uniform bb',128,cTemp);
scene = sceneAdjustLuminance(scene,lum);
scene = sceneSet(scene,'fov',10);
oi = oiCompute(oi,scene,'pad value','mean');
oi = oiCrop(oi,'border');
% oiShowImage(oi);

cm.compute(oi);

%% Calculate the absorptions

lCones = conerectGet(cm,'absorptions','l cones');
mCones = conerectGet(cm,'absorptions','m cones');
sCones = conerectGet(cm,'absorptions','s cones');
fprintf('At (%.1f,%.1f) deg position, absorptions per candela: %.1f,%.1f,%.1f \n',cm.center, [mean(lCones), mean(mCones), mean(sCones)])

% Rough summary across cone
fprintf('Absorptions per candela/m^2  %.1f\n',[mean(lCones),mean(mCones),mean(sCones)]/lum)

%% END