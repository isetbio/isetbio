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
    % Default TreeShrew model
    oi = oiTreeShrewCreate();
%}
%{
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
    
    if checkfields(oi.optics, 'transmittance')
        oi.optics = rmfield(oi.optics, 'transmittance');
    end  
    
    % TreeShrew lens absorption.
    theTreeShrewLens = treeShrewLens(oiGet(oi, 'optics wave'));
    
    % Update the oi.lens
    oi = oiSet(oi, 'lens', theTreeShrewLens);
    
    % Set the same lens in the optics structure too
    oi.optics.lens = theTreeShrewLens;
end

% Method to generate a Lens object with treeshrew - based optical density
function theTreeShrewLens = treeShrewLens(targetWavelenth)
    % TreeShrew lens absorption. Start with the human lens.
    theTreeShrewLens = Lens();
    
    % Load TreeShrew lens unit-density
    load('treeshrewLensAbsorbance.mat', 'wavelength', 'data');
    
    % Interpolate to optics wavelength support
    unitDensity = interp1(wavelength,data,targetWavelenth, 'pchip');
    
    % Update the lens
    set(theTreeShrewLens,'wave', targetWavelenth);
    set(theTreeShrewLens,'unitDensity',unitDensity);
end