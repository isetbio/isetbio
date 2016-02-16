function receptiveFieldDiameterParasol2STD = receptiveFieldDiameterFromTEE(tee)
%% receptiveFieldDiameter2STD converts the TEE to RF diameter 
% 
%           receptiveFieldDiameterParasol2STD = receptiveFieldDiameterFromTEE(temporalEquivEcc);
%               [only called internally from @rgcMosaic/initialize.m]
% 
% The RF diameter is in um at the specified TEE. The RF
% diameter can then be used to bulid RFs for the rgc object in
% buildSpatialRFArray.m.
% 
% See Chichilnisky, E. J., and Rachel S. Kalmar. "Functional asymmetries 
% in ON and OFF ganglion cells of primate retina." The Journal of 
% Neuroscience 22.7 (2002), Fig. 5, pg. 2741.
% 
% Example:
%   receptiveFieldDiameterParasol2STD = receptiveFieldDiameterFromTEE(innerRetina.temporalEquivEcc);
% 
% See also: @rgcMosaic/initialize.m

%% Define numerical parameters from published data
% Fig. 5, Chichilnisky & Kalmar 2002 shows dendritic field diameter as a
% function of TEE; DF diameter = 1.57*(RF diameter)
scaleFactor = 1.57; %

ecc = [0.5 10]; dia2STD = [25  275]/scaleFactor;

%% Estimate linear fit to points in Fig. 5 of C&K 2002
% DF diameter = m*TEE + yint
m = (dia2STD(2)-dia2STD(1))/(ecc(2)-ecc(1));
yint = dia2STD(1) - m*ecc(1);

%% Get our RF diameter for a particular TEE
receptiveFieldDiameterParasol2STD = m*tee + yint;