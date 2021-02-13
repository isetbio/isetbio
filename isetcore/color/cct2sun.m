function spd = cct2sun(wave, cct, units)
% Correlated color temperature to daylight SPD at specified wavelengths
%
% Syntax:
%   spd = cct2sun( [wave], cct, [units])
%
% Description:
%    Determines daylight/sun spectral power distribution based on
%    correlated color temperature.
%
% Inputs:
%    wave  - (Optional) Vector. Wavelength vector of SPD. Default 400:700.
%    cct   - Numeric. Correlated color temperatures.
%            [Note: XXX - Can be a vector, but should it?]
%    units - (Optional) String. A string that specifies units of returned
%            spd (spectral power distribution). Default 'energy'. Options:
%       'energy': Default. Also 'watts', but please use 'energy'
%       'photons': Also 'quanta', but please use 'photons'
%
% Outputs:
%    spd   - Vector. The daylight or sun SPD in units of (relative) energy
%            or photons.
%
% Optional key/value pairs:
%    None.
%
% References:
%    http://en.wikipedia.org/wiki/Standard_illuminant
%    Judd, Macadam, Wyszecki
%       http://www.opticsinfobase.org/abstract.cfm?URI=josa-54-8-1031
%
% See Also:
%   daylight
%

% History:
%    08/14/00  xxx  Last Updated
%    10/30/17  jnm  Comments & formatting
%    11/17/17  jnm  Note, final example & formatting
%    12/12/17  bw   Checked and removed daylight basis comment, which
%                   was fixed
%    07/11/19  JNM  Formatting update

% Examples:
%{
    wave = 400:5:700;
    cct = 4000;
    spd = cct2sun(wave, cct);
    plot(wave, spd);
%}
%{
    wave = 400:2:700;
    cct = 6500;
    spd = cct2sun(wave, cct, 'photons');
    plot(wave, spd);
%}
%{
    wave = 400:2:700;
    cct = 6500;
    spd = cct2sun(wave, cct, 'energy');
    plot(wave, spd);
%}
%{
    w = 400:700;
    spd = cct2sun(w, [4000 6500], 'photons');
    plot(w, spd)
%}
%{
    spd = cct2sun([],[4000:5000]);
    plot(400:700, spd)
%}

if notDefined('wave'), wave = 400:700; end
if notDefined('cct'), error('Correlated color temperature required'); end
if notDefined('units'), units = 'energy'; end

% Calculate the xy chromaticity coordinates.
mask = 1 .* (cct >= 4000 & cct < 7000 ) + 2 .* (cct >= 7000 & cct < 30000);

ind = find(mask == 0, 1);
if ~isempty(ind)
   error('At least one CCT is outside the acceptable range [4000-30000]');
end

% Look this up and put in a reference to the appropriate W&S pages.
xdt = zeros(2, size(cct, 2));
xdt(1, :) = -4.6070e9 ./ cct .^ 3 + 2.9678e6 ./ cct .^ 2 + 0.09911e3 ./ ...
    cct + 0.244063;
xdt(2, :) = -2.0064e9 ./ cct .^ 3 + 1.9018e6 ./ cct .^ 2 + 0.24748e3 ./ ...
    cct + 0.237040;

% Explain the mask terms.
xd = (mask == 1) .* xdt(1, :) + (mask == 2) .* xdt(2, :);
yd = -3.000 * xd .^ 2 + 2.870 * xd - 0.275;

% Calculate the CIE SUN weights that will be applied to the CIE sunlight
% basis functions in CIESUN.
M  = zeros(2, size(cct, 2));
M(1, :) = (-1.3515 - 1.7703 * xd + 5.9114 * yd) ./ ...
   (0.0241 + 0.2562 * xd - 0.7341 * yd);
M(2, :) = (0.0300 - 31.4424 * xd + 30.0717 * yd) ./ ...
   (0.0241 + 0.2562 * xd - 0.7341 * yd);

% Calculate the final daylight SPD.
dayBasis = ieReadSpectra('cieDaylightBasis', wave);
spd = dayBasis(:, 2:3) * M + repmat(dayBasis(:, 1), [1 size(cct, 2)]);

% Flip to photons/quanta if needed. Energy/watts would be the alternative.
switch lower(units)
    case {'quanta', 'photons'}
        spd = Energy2Quanta(wave, spd);
    case {'energy', 'watts'}
        % Don't need to do anything here.
    otherwise
        error('Unknown units specified');
end

end