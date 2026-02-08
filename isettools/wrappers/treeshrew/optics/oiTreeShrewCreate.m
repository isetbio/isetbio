function oi = oiTreeShrewCreate(varargin)
% Create an optical image structure specific to the TreeShrew
%
% Syntax:
%   oi = OITREESHREWCREATE(oiType, [varargin]) 
%
% See Also:
%    opticsTreeShrewCreate
%
% History:
%    11/23/18  NPC  ISETBIO TEAM, 2018
%
% Examples:
%{
    % ETTBSkip
    % This is currently broken because we never updated after the big 2024 
    % ISETBio/ISETCam consilidation.

    % Default TreeShrew model
    oi = oiTreeShrewCreate();
%}
%{
    % ETTBSkip
    % This is currently broken because we never updated after the big 2024 
    % ISETBio/ISETCam consilidation.

    % Custom tree shew optics where we specify the sigma of the PSF
    oi = oiTreeShrewCreate(...
        'inFocusPSFsigmaMicrons', 40 ... % 40 microns
    );
%}

    oi.type = 'opticalimage';
    oi = oiSet(oi, 'optics', opticsTreeShrewCreate(varargin{:}));
    oi = oiSet(oi, 'name', 'treeshrew');
    
    oi = oiSet(oi, 'bit depth', 32);
    oi = oiSet(oi, 'diffuser method', 'skip');
    oi = oiSet(oi, 'consistency', 1);
    
    % TreeShrew lens absorption.  For ISETCam integration, need to 
    % remove transmittance field first, which we do by brute force.
    if checkfields(oi.optics, 'transmittance')
        oi.optics = rmfield(oi.optics, 'transmittance');
    end
    theTreeShrewLens = lensTreeShrewCreate('wave', oiGet(oi, 'optics wave'));
    
    % Update the oi.lens
    oi = oiSet(oi, 'optics lens', theTreeShrewLens);
    
    % Set the same lens in the optics structure too
    oi.optics.lens = theTreeShrewLens;
end