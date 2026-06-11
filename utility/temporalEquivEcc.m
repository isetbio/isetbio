function tee = temporalEquivEcc(center, varargin)
% (IN PROG.) Equivalent temporal eccentricity acc. to Kalmar/Chichilnisky
%
% Syntax:
%   tee = temporalEquivEcc(center)
%
% Description:
%    Calculate the temporal equivalent eccentricity using formula from
%    Chichilnisky, E. J., and Rachel S. Kalmar. "Functional
%    asymmetries in ON and OFF ganglion cells of primate retina." The
%    Journal of Neuroscience 22.7 (2002).
%
%    This function compensates for the observation that
%    "contour lines of constant RGC density are circular in the temporal
%    half but elliptical in the nasal half of the retina" (Chichilnisky &
%    Kalmar, pg. 2738, 2002).
%
%    The TEE can be used, for example, to calculate the diameter of the
%    receptive field size at a retinal location.
%      TEE = sqrt((0.61 * X ^ 2) + Y ^ 2) (corrected from the publication)
%
%    Examples contained in code.
%
% Inputs:
%    center  - position on retina in mm, fovea is (0, 0)
%
% Outputs:
%    tee     - Calculated temporal equivalent eccentricity
%
% Optional key/value pairs:
%    'eyeSide' - Which eye? (Default 'left'). Not yet implemented in code.
%
% Notes:
%    * [Note: DHB - This routine is not really done.  No implementation of
%      which eye yet, and what the coordinate system is, is not specified.
%      And, the current code doesn't take the quadrant of the retina into
%      account. Maybe the commented out code in the note from BAW  below is
%      right.  All needs to be fixed up before we use this for anything
%      serious. It is only called from the rgc method initSpace.  I put in
%      an error message on execution, as this seemed really broken.]
%    * [Note: DHB - Key 'eyeSide' is not well named, 'whichEye' would be
%      better.]
%    * [Note: BAW - Below is the from JRG. Not sure it is all the same.
%      In general, I am not sure we should use this for human work.]

%     % Convert angle in degrees to radians
%     theta = (pi / 180) * theta;
%
%     %% Apply formula for TEE
%     if ((theta > (pi / 2) && theta < (3 * pi / 2)) && eyeSide == 1) ...
%             || ((theta < (pi / 2) || theta > (3 * pi / 2)) ...
%             && eyeSide == 0)
%
%         [xrad, yrad] = pol2cart(theta, rho);
%
%         % Apparently some correction for other quadrants. Formula below:
%         aspectRatio = 0.61;
%         TEE = sqrt((xrad * aspectRatio) .^ 2 + yrad .^ 2);
%
%     else
%         % Temporal side. No need for a correction.
%         TEE = rho;
%     end

% History:
%    xx/xx/17   BW  ISETBIO Team, 2017
%    11/20/17  jnm  Comments & formatting
%    01/24/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    % ETTBSkip - This example known not to work.  Remove line when fixed.
    temporalEquivEcc([1, -1],'eyeSide', 'left')
    temporalEquivEcc([-1, 1],'eyeSide', 'left')
%}

error('This routine is too unfinished to use. See Notes in source.');


%%
p = inputParser;
p.addRequired('center', @isvector);
p.addParameter('eyeSide', 'left', @ischar);
p.parse(center, varargin{:});


% Not handled yet. See above.
% eyeSide = p.Results.eyeSide;

tee = sqrt((0.61 * center(1) ^ 2) + center(2) ^ 2);

end