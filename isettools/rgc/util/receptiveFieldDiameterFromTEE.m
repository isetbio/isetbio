function receptiveFieldDiameterParasol2STD = receptiveFieldDiameterFromTEE(tee)
%% receptiveFieldDiameter2STD converts the TEE to RF diameter 
% 
%           receptiveFieldDiameterParasol2STD = receptiveFieldDiameterFromTEE(temporalEquivEcc);
%               [only called internall from @rgcMosaic/initialize.m]
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

scaleFactor = 1.57; % DF diameter = 1.57*(DF diameter)

ecc = [0.5 10]; dia2STD = [25  275]/scaleFactor;
m = (dia2STD(2)-dia2STD(1))/(ecc(2)-ecc(1));
yint = dia2STD(1) - m*ecc(1);

receptiveFieldDiameterParasol2STD = m*tee + yint;