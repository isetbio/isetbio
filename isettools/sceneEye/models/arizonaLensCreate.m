function [az, txt] = arizonaLensCreate(accommodation)
% Return the parameters needed to write the Arizona eye model
%
% Synopsis
%  [az, txt] = arizonaLensCreate(accommodation)
%
% Inputs
%  accommodation:  In diopters (1/m)
%
% Optional key/value p;airs
%   N/A
%
% Returns
%   az:   Struct with the relevant rows of the lens matrix.
%   txt:  Text description of the columns
%
% Description
%  https://photonengr.com/wp-content/uploads/kbasefiles/ArizonaEyeModel.pdf
%
% See also
%  arizonaWrite, navarroWrite, legrandWrite
% 

% Examples:
%{
   [az, columnDescription]  = arizonaLensCreate(0);
%}

%% Set the basic model parameters

txt = 'radiusX, radiusY, thickness, materialIndex, semiDiameter, conicConstantX, conicConstantY';

% The following equations define dependences based on accommodation
% (diopters).
%
% They are taken from:
% 
%   https://photonengr.com/wp-content/uploads/kbasefiles/ArizonaEyeModel.pdf
%
anteriorRadius   = 12 - (0.4 * accommodation);
posteriorRadius  = -5.224557 + (0.2 * accommodation); 
aqueousThickness = 2.97 - (0.04 * accommodation);
lensThickness    = 3.767 + (0.04 * accommodation);
anteriorAsph     = -7.518749 + (1.285720 * accommodation);
posteriorAsph    = -1.353971 - (0.431762 * accommodation);

%% 
% Columns are: 
%          radiusX, radiusY, thickness, materialIndex, semiDiameter, conicConstantX, conicConstantY
corneaA = [-7.8,    -7.8,      0.55,   1,     4.820,    -0.25,     -0.25];
corneaP = [-6.5,    -6.5,      2.97,   2,     4.341,     -0.25,    -0.25];
pupil   = [0,        0,        0,      2,     2,          0,        0];
lensA   = [-12,     -12,       3.767,  3,     3.750,    -7.518749, -7.518749];
lensP   = [5.224557, 5.224557, 0,      4,     3.750,    -1.353971, -1.353971];

% flip these because of sign conventions
lensA(1:2) = [anteriorRadius anteriorRadius]   * -1;
lensP(1:2) = [posteriorRadius posteriorRadius] * -1;

% Accommodation adjustments
corneaP(3) = aqueousThickness;
lensA(3)   = lensThickness;
lensA(6:7) = [anteriorAsph anteriorAsph];
lensP(6:7) = [posteriorAsph posteriorAsph];

% Build the struct for simpler return
az.corneaA = corneaA;
az.corneaP = corneaP;
az.pupil   = pupil;
az.lensA   = lensA;
az.lensP   = lensP;

end