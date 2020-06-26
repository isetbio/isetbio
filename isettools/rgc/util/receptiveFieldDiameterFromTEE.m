function receptiveFieldDiameterParasol2STD = ...
    receptiveFieldDiameterFromTEE(tee)
% Convert the TEE (mm) to RF diameter (um)
%
% Syntax:
%   receptiveFieldDiameterParasol2STD = ...
%       receptiveFieldDiameterFromTEE(tee);
%
% Description:
%    I think the TEE is in mm on the retinal surface.  The function to
%    compute it is
%
%          TEE = retinalLocationToTEE(theta, rho, eyeSide)
%
%    The RF diameter is in um at the specified TEE. The RF diameter can
%    then be used to bulid RFs for the rgc object in buildSpatialRFArray.m.
% 
%    See Chichilnisky, E. J., and Rachel S. Kalmar. "Functional asymmetries
%    in ON and OFF ganglion cells of primate retina." The Journal of
%    Neuroscience 22.7 (2002), Fig. 5, pg. 2741.
%
%    This function contains examples of usage inline. To access them, type
%    'edit receptiveFieldDiameterFromTEE.m' into the Command Window.
%
% Inputs:
%    tee - Numeric. The temporal Equivalent Eccentricity.
%
% Outputs:
%    receptiveFieldDiameterParasol2STD
%        - Numeric. The receptive field diameter.
%
% Optional key/value pairs:
%    None.
%

% History:
%    XX/XX/17  BW   ISETBIO Team, 2017
%    06/05/19  JNM  Documentation pass

% Examples:
%{
    % ETTBSkip - skipping broken examples
    %
    % This needs code to define rgcLayer before it could
    % possibly work.
    TEE = retinalLocationToTEE(rgcLayer.center, eyeSide)
    receptiveFieldDiameterParasol2STD = receptiveFieldDiameterFromTEE(TEE);
%}

%% Define numerical parameters from published data
% Fig. 5, Chichilnisky & Kalmar 2002 shows dendritic field diameter as a
% function of TEE; DF diameter = 1.57 * (RF diameter)
scaleFactor = 1.57; %
ecc = [0.5 10];
dia2STD = [25  275] / scaleFactor;

%% Estimate linear fit to points in Fig. 5 of C&K 2002
% DF diameter = m * TEE + yint
m = (dia2STD(2) - dia2STD(1)) / (ecc(2) - ecc(1));
yint = dia2STD(1) - m * ecc(1);

%% Get our RF diameter for a particular TEE
receptiveFieldDiameterParasol2STD = m * tee + yint;

end
