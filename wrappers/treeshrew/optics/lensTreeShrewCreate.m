function theLens = lensTreeShrewCreate(varargin)
% Create a Lens object with treeshrew - based optical density
%
% Syntax:
%   lens = LENSTREESHREWCREATE([varargin]) 
%
% Inputs: 
%           None
%
%
% Optional key/value pairs:
%    wave                 - Vector. Wavelengths. Default [], retrieve from
%                                   absorbance file
%    lensAbsorbanceFile   - Char, Filename containing 'wavelength' and
%                                   'data' fields specifying the lens
%                                   absorbance spectrum
%
% Outputs: 
%           lens - Struct. The created lens structure.
%
% See Also:
%    opticsTreeShrewCreate
%
% History:
%    11/23/18  NPC  ISETBIO TEAM, 2018
%
% Examples:
%{
    % Default TreeShrew lens model
    lens = lensTreeShrewCreate();
%}

    %% parse input
    p = inputParser;
    p.addParameter('wave', [], @isnumeric);
    p.addParameter('lensAbsorbanceFile', 'treeshrewLensAbsorbance.mat', @ischar);
    p.parse(varargin{:});
    
    lensAbsorbanceFile = p.Results.lensAbsorbanceFile;
    targetWavelenth = p.Results.wave;
    
    % TreeShrew lens absorption. Start with the human lens.
    theLens = Lens();
    
    % Load TreeShrew lens unit-density
    load(lensAbsorbanceFile, 'wavelength', 'data');
    
    if (isempty(targetWavelenth))
        targetWavelenth = wavelength;
        unitDensity = data;
    else
        % Interpolate to optics wavelength support
        unitDensity = interp1(wavelength,data,targetWavelenth, 'pchip');
    end
    
    
    % Update the lens
    set(theLens,'wave', targetWavelenth);
    set(theLens,'unitDensity',unitDensity);
end