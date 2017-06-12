function ecc = temporalEquivEcc(center,varargin)
% Equivalent temporal eccentricity according to Kalmar/Chichilnisky
%
% Calculate spatial RF diameter for ON Parasol cell at a particular TEE
% See Chichilnisky, E. J., and Rachel S. Kalmar. "Functional asymmetries
% in ON and OFF ganglion cells of primate retina." The Journal of
% Neuroscience 22.7 (2002), Fig. 5, pg. 2741. 2STD fit in micrometers.
% This function compensates for the observation that "contour lines of
% constant RGC density are circular in the temporal half but elliptical in
% the nasal half of the retina" (Chichilnisky & Kalmar, pg. 2738, 2002).
%
% The TEE is used to calculate the diameter of the receptive field size at
% that location, given by Fig. 5 Chichilnisky & Kalmar (2002).
%
%     TEE = sqrt((0.61*X^2)+Y^2) (corrected from the publication)
%
% BW ISETBIO Team, 2017
%
% This is the code from JRG.  We need to deal with parameters
%
%% Get eye side
% if strcmp(eyeSide,'right'),        eyeSide = 0;
% elseif strcmp(eyeSide,'left'),     eyeSide = 1;
% else
%     error('Incorrect eye side specification %s\n',eyeSide);
% end
% 
% % Convert angle in degrees to radians
% theta = (pi/180)*theta;
% 
% %% Apply formula for TEE
% if ( (theta > (pi/2) && theta < (3*pi/2)) && eyeSide==1 ) || ...
%         ( (theta < (pi/2) || theta > (3*pi/2)) && eyeSide==0 )
%     
%     [xrad, yrad] = pol2cart(theta, rho);
%     
%     % Apparently some correction for other quadrants.  Here is the formula.
%     aspectRatio = 0.61;
%     TEE = sqrt((xrad*aspectRatio).^2 + yrad.^2);
%     
% else
%     % Temporal side.  No need for a correction.
%     TEE = rho;
% end

p = inputParser;
p.addRequired('center',@isvector);
p.addParameter('eyeSide','left',@ischar);
p.parse(center,varargin{:});

% Not handled yet.  See above.
% eyeSide = p.Results.eyeSide;

ecc = sqrt((0.61*center(1)^2)+center(2)^2);

end