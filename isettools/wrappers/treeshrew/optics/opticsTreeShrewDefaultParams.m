function defaultParams = opticsTreeShrewDefaultParams()
% Generate the default parameter values for the tree-shrew eye
%
% Syntax:
%   opticsParams = OPTICSTREESHREWDEFAULTPARAMS(varargin)
%
% Default params:
% The peak spatial frequency sensitivity is around 0.5 cyc/deg, 
% with a high-frequency cutoff of 2 c/deg (see "Spatial contrast sensitivity of the tree shrew". Petry HM, Fox R, Casagrande VA.)
% 
% From "Normal development of refractive state and ocular component dimensions 
% in the tree shrew (Tupaia belangeri)", by Thomas T. Norton, Neville A. McBrien, 1992
% anteriorFocalLengthMM = 4.35;
% posteriorNodalDistanceMM = posteriorFocalPoint - posteriorNodalPoint = 7.84-3.49 = 4.35;
%
defaultParams = struct(...
    'opticsType', 'gaussian psf', ...
    'inFocusPSFsigmaMicrons', 7.0, ...
    'focalLengthMM', 4.35, ...
    'pupilDiameterMM', 3.0 ...
    );

end
