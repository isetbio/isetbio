function set(lens, param, val, varargin)
% Set parameters of the lens structure
%
%     Lens.set(param,varargin)
%
% See Lens.get for notes about the properties in general and the
% formulae relating absorbance and absorptance and transmittance.
%
% Inputs:
%   lens   - Lens class object
%   param  - string, paramter to be set
%   val    - desired value
%
% param can be chosen from
%   name          - object name
%   wave          - wavelength samples
%   unitDensity   - absorbance data
%   density       - single value
%
% Notes:
%   Absorbance spectra are normalized to a peak value of 1.
%   Absorptance spectra are the proportion of quanta actually absorbed.
%   Equation: absorptanceSpectra = 1 - 10.^(-OD * absorbanceSpectra)
%
% HJ/BW, ISETBIO TEAM, 2016

% parse input
p = inputParser; p.KeepUnmatched = true;
p.addRequired('param', @isstr);
p.addRequired('val');
p.parse(param, val, varargin{:});

% set parameter value
switch ieParamFormat(param)
    case 'name'
        lens.name = val;
    case {'wave', 'wavelength'}
        % Sampling wavelength
        lens.wave = val(:);
    case {'absorbance','unitdensity'}
            assert(length(val) == length(lens.wave), ...
                'Val should have same length as lens wavelength');
            lens.unitDensity = val;
    case 'density'
        assert(isscalar(val), 'val should be scalar');
        lens.density = val;      
    otherwise
        error('Unknown parameter %s\n',param);
end

end