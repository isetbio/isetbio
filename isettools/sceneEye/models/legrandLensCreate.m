function lg = legrandLensCreate
% Returns parameters of the LeGrand eye model
%
% Synopsis:
%  lg = legrandWrite()
%
% Inputs
%   N/A
%
% Optional key/value pairs
%   N/A
%
% Outputs
%   lg  - Struct with the parameters for the different radii, curvatures
%         and such of the model
%   ior - Matrix with the indices of refraction
%
% Description
%   The LeGrand eye model does not have a set for accommodation or a model
%   for how to change the accommodation.  The Arizona eye and the Navarro
%   eye do have set functions (see below)
%
% See also
%   legrandWrite, navarroWrite,setNavarroAccommodation, setArizonaAccommodation
%

%%  Source
%
% Atchison, David A., and George Smith. "Chromatic dispersions of the
% ocular media of human eyes." JOSA A 22.1 (2005): 29-37.

%% Here are the geometrical parameters

%             radiusX, radiusY, thickness, materialIndex, semiDiameter, conicConstantX, and conicConstantY
lg.corneaA = [-7.8,     -7.8, 0.55, 1, 4.820, 0, 0];
lg.corneaP = [-6.5,     -6.5, 3.05, 2, 4.341, 0, 0];
lg.pupil   = [0,           0,  0,   2, 2,     0, 0];
lg.lensA   = [-10.2,   -10.2,  4,   3, 3.750, 0, 0];
lg.lensP   = [6,          6,   0,   4, 3.750, 0, 0];

end
