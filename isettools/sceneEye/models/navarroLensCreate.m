function [na,txt] = navarroLensCreate(accommodation)
% Create a data struct for the Navarro eye model
%
% Syntax:
%   [na,txt] = navarroLensCreate(accommodation)
%
% Description:
%    Create a Navarro eye model data structure.
%
% Inputs:
%    accommodation: The accommodation in diopters
%
% Optional key/val pairs
%    N/A
%
% Outputs:
%    None.
%
% See also
%   setNavarroAccommodation

% Examples:
%{
 [na,txt] = navarroLensCreate(0);
%}

%% Check the variable
if notDefined('accommodation'), error('Accommodation must be defined'); end

txt = 'radiusX, radiusY, thickness, materialIndex, semiDiameter, conicConstantX, and conicConstantY';

%% These equations are from Table 4 in Navarro's paper.
anteriorRadius   = 10.2 - 1.75 * log(accommodation + 1);
posteriorRadius  = -6 + 0.2294 * log(accommodation + 1);
aqueousThickness = 3.05 - 0.05 * log(accommodation + 1);
lensThickness    = 4 + 0.1 * log(accommodation + 1);
anteriorAsph     = -3.1316 - 0.34 * log(accommodation + 1);
posteriorAsph    = -1 - 0.125 * log(accommodation + 1);

% Columns are: 
%  radiusX, radiusY, thickness, materialIndex, semiDiameter, conicConstantX, and conicConstantY
corneaA = [-7.72, -7.72, 0.55, 1, 4.820, -0.26, -0.26];
corneaP = [-6.5, -6.5, 3.050, 2, 4.341, 0, 0];
pupil   = [0, 0, 0, 2, 2, 0, 0];
lensA   = [-10.2, -10.2, 4, 3, 3.750, -3.132, -3.132];
lensP   = [6, 6, 0, 4, 3.750, -1, -1];

% flip these because of sign conventions
lensA(1:2) = [anteriorRadius anteriorRadius] * -1;
lensP(1:2) = [posteriorRadius posteriorRadius] * -1;

corneaP(3) = aqueousThickness;
lensA(3)   = lensThickness;
lensA(6:7) = [anteriorAsph anteriorAsph];
lensP(6:7) = [posteriorAsph posteriorAsph];

%% Build the struct for the return
na.corneaA = corneaA;
na.corneaP = corneaP;
na.pupil   = pupil;
na.lensA   = lensA;
na.lensP   = lensP;

end
