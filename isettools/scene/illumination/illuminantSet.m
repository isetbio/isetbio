function il = illuminantSet(il, param, val, varargin)
% Set parameter value for illuminant structure
%
% Syntax:
%   il = illuminantSet(il, param, val, varargin)
%
% Description:
%    The illuminant structure has various formats. It can be a simple
%    vector, which defines an SPD for the entire image. It can also be an
%    RGB format spectral data set that has a separate illuminant SPD at
%    every point.
%
%    We are in the middle of extending this spatial-spetral format and
%    writing tutorials.
%
% Inputs:
%    il       - The illuminant to modify
%    param    - The parameter to change
%    val      - The value to attribute to param
%    varargin - (Optional) Addition variables as required. Default is []
%
% Outputs:
%    il       - The modified illuminant
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note: JNM - The second example works, but I'm not sure that I've
%      provided realistic values for photons or energy.]
%
% See Also:
%    illuminantCreate, illuminantGet
%

% History:
%    xx/xx/12       (c) ImagEval Consulting, LLC, 2012.
%    01/29/18  jnm  Formatting, fix second example.
%    03/26/18  dhb  Fix type (can't write *e-3, need to write *1e-3).

% Examples:
%{
    % Create the illuminant structure, with no data
	il = illuminantCreate;
    il = illuminantSet(il, 'name', 'outdoor');
%}
%{
    il = illuminantCreate;
    photons = 256;
    e = zeros(1, length(illuminantGet(il, 'wave')));
    il = illuminantSet(il, 'photons', photons);
    % or the following will convert to energy for you
    il = illuminantSet(il, 'energy', e);
%}

%%
if notDefined('il'),   error('illuminant structure required'); end
if notDefined('param'), error('param required'); end
if ~exist('val', 'var'), error('val is required'); end

%%
param = ieParamFormat(param);
switch param
    case 'name'
        il.name = val;
    case 'type'
        if ~strcmpi(val, 'illuminant')
            error('Type must be illuminant');
        end
        il.type = val;
    case 'photons'
        % il = illuminantSet(il, 'photons', data);
        il.data.photons = single(val);
    case 'energy'
        % User sent in energy. We convert to photons and set.
        % We need to handle the spatial spectral case properly.
        % See s_illuminantSpace
        wave = illuminantGet(il, 'wave');
        if ndims(val) > 2 %#ok<ISMAT>
            [val, r, c] = RGB2XWFormat(val);
            val = Energy2Quanta(wave, val')';
            val = XW2RGBFormat(val, r, c);
            il = illuminantSet(il, 'photons', val);
        else
            % For set of a vector to be a column vector
            il = illuminantSet(il, 'photons', Energy2Quanta(wave, val(:)));
        end
    case {'wave', 'wavelength'}
        % il = illuminantSet(il, 'wave', wave)
        % Need to interpolate data sets and reset when wave is adjusted.
        oldW = illuminantGet(il, 'wave');
        newW = val(:);
        il.spectrum.wave = newW;

        % Now decide what to do with photons
        p = illuminantGet(il, 'photons');
        if ~isempty(p)
            % If p has the same length as newW, let's assume it was already
            % changed. Otherwise, if it has the length of oldW, we should
            % try to interpolate it.
            if length(p) == length(newW)
                % Sample length of photons already equal to newW. No 
                % problems here.
            elseif length(p) == length(oldW)
                % Adjust the sampling.
                newP = interp1(oldW, p, newW, 'linear', min(p(:) * 1e-3)');
                il = illuminantSet(il, 'photons', newP);
            else 
                error(['Photons and wavelength sample points not '
                    'interpretable']);
            end
            % vcNewGraphWin; plot(newW, newP);
        end
    case 'comment'
        il.comment = val;
    otherwise
        error('Unknown illuminant parameter %s\n', param)
end

end
