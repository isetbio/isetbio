function XYZ = ieXYZFromEnergy(energy, wave)
% CIE XYZ values from spectral radiance(w/nm/sr/m2) or irradiance(w/nm/m2)
%
% Syntax:
%   XYZ = ieXYZFromEnergy(energy, wave)
%
% Description:
%    Calculate the XYZ values of the spectral radiance or irradiance
%    functions in the variable energy. The input format of energy can be
%    either XW (space-wavelength) or RGB. The wavelength samples of energy
%    are stored in the variable wave.
%
%    Notice, that XW is AN UNUSUAL FORMAT for energy. Often, we put the
%    SPDs into the columns of the matrix. But in the XW format, the SPDs
%    are in the rows. Sorry.
%
%    The returned values, XYZ, are X, Y, Z in the columns of the matrix. Each
%    row of energy has a corresponding XYZ value in the corresponding row
%    of XYZ. This is what we call XW format.
%
%    The units of Y are candelas/meter-squared if energy is radiance and
%    lux if energy is irradiance.
%
%    The units of radiance are watts/[sr-m2-nm], while the units of
%    irradiance are Watts/[m2-nm].
%
% Inputs:
%    energy - Energy in XW (space-wavelength) or RGB formats
%    wave   - Wavelength samples of energy
%
% Outputs:
%    xyz    - X, Y, and Z are columns in matrix format. Y is in candelas
%             per meter squared (radiance), or lux (irradiance). Units for
%             radiance are watts per the combination of steradian times
%             meters squared times nanometer (w/[sr*m^2*nm]), and for
%             irradiance are watts per the combination of meters squared
%             times nanometers (w/[m^2*nm])
%
% Notes:
%    * We could consider putting the return into RGB format if it is sent
%      in that way.
%    * XXX - Returning in RGB forma is new. I tested it in the scielab
%      branch with v_ISET and some other calls. But it might cause
%      something to break somewhere. Stay alert!
%
% See Also:
%   ieXYZFromPhotons()
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC.
%    10/27/17  jnm  Comments & formatting

% Examples:
%{
   wave = 400:10:700;
   tmp = load('CRT-Dell');
   dsp = tmp.d;
   energy = displayGet(dsp, 'spd', wave);
   energy = energy';
   displayXYZ = ieXYZFromEnergy(energy, wave)

   patchSize = 1;
   macbethChart = sceneCreate('macbeth', patchSize); 
   p = sceneGet(macbethChart, 'photons');
   wave = sceneGet(macbethChart, 'wave');
   e = Quanta2Energy(wave, p);
   XYZ = ieXYZFromEnergy(e, wave);  
%}

% Force data into XW format.
if ndims(energy) == 3
    if length(wave) ~= size(energy, 3)
        error('Bad format for input variable energy.');
    end
end

% Returning in RGB forma is new. I tested it in the scielab branch with
% v_ISET and some other calls. But it might cause something to break
% somewhere. Stay alert!
iFormat = vcGetImageFormat(energy, wave);
switch iFormat
    case 'RGB'
        % [rows, cols, w] = size(data);
        [xwData, r, c] = RGB2XWFormat(energy);
        % disp('RGB return')
    otherwise
        % XW format
        xwData = energy;
end

% xwData = ieConvert2XW(energy, wave);
if size(xwData, 2) ~= length(wave)
    error('Problem converting input variable energy into XW format.');
end

% The spectra of the energy points are in the rows of xwData. We load the
% XYZ color matching functions into the columns of S.
S = ieReadSpectra('XYZ', wave);
if numel(wave) > 1
    dWave = wave(2) - wave(1);
else
    dWave = 10;
    disp('10 nm band assumed');
end

% The return value has three columns, [X, Y, Z].
XYZ = 683 * dWave * (xwData*S);

% If it was sent in RGB, return it in RGB
switch iFormat
    case 'RGB'
        XYZ = XW2RGBFormat(XYZ, r, c);
    otherwise
        % XW format, do nothing
end

end
