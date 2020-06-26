function [illuminance, meanIlluminance, meanCompIlluminance] = ...
    oiCalculateIlluminance(oi)
% Calculate illuminance (lux) of optical image spectral irradiance.
%
% Syntax:
%   [illuminance, meanIlluminance, meanCompIlluminance] = ...
%       oiCalculateIlluminance(opticalImage)
%
% Description:
%    The optical image spectral irradiance data are converted into
%    illuminance (Lux) using the CIE formula.
%
%    Suppose the spectral irradiance is irradianceE (watts/m2) and sampled
%    at various wavelength values (nm); vLambda is the photopic sensitivity
%    function sampled at the same set of wavelengths; suppose the
%    wavelength spacing is binwidth (nm). Then the formula for illuminance
%    in units of lux is
%
%       illuminance = (683 * binwidth) * irradianceE * vLambda;
%
%    The mean illuminance can also be computed and returned. The
%    complementary illuminance is computed using (1 - v(lambda)) rather
%    than v(lambda). This is calculated for certain infrared applications.
%
%    There are examples contained in the code below. To access them, type
%    'edit oiCalculateIlluminance.m' in to the Command Window.
%
% Inputs:
%    oi                  - Struct. An optical image structure.
%
% Outputs:
%    illuminance         - Matrix. The illuminance of the optical image's
%                          irradiance in lux.
%    meanIlluminance     - Numeric. The mean illuminance of the OI.
%    meanCompIlluminance - The complementary illuminance
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * TODO: There is something odd with the compIlluminance calculation
%      below. Need to discuss with MP.
%    * TODO: We need to plug in a function below that describes a new
%      metric for quantifying the effect of IR energy. (m.p. 12/16/07)
%    * TODO: Determine if we want examples that have illuminances that are
%      not blank/empty?
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    03/02/18  jnm  Formatting
%    06/26/19  JNM  Documentation update, added third example, uncommented
%                   newPeak to fix failing script if meanCompIlluminance is
%                   requested. Not sure if 750 value is the correct peak to
%                   specify, so I'll leave that up to brighter minds.

% Examples:
%{
    oi = oiCreate;
    illuminance = oiCalculateIlluminance(oi);
    oi = oiSet(oi, 'illuminance', illuminance);
%}
%{
   % Complementary illuminance mean
   oi = oiCreate;
   [I, meanI, meanC] =  oiCalculateIlluminance(oi);
%}
%{
    scene = sceneCreate('uniform');
    scene = sceneSet(scene, 'fov', 15);  % Reasonably large
    scene = sceneAdjustLuminance(scene, 10 ^ -10);

    oi = oiCreate;
    % No lens shading
    optics = oiGet(oi, 'optics');
    optics = opticsSet(optics, 'cos4th', 'off');
    oi = oiSet(oi, 'optics', optics);
    oi = oiCompute(oi, scene);

   [I, meanI] =  oiCalculateIlluminance(oi);
%}

if notDefined('oi'), error('Optical image required.'); end

wave = oiGet(oi, 'wave');
binWidth = oiGet(oi, 'binWidth');
sz = oiGet(oi, 'size');

% Infrared (complementary) illuminance
meanCompIlluminance = 0;

% Read the V-lambda data at the relevant wavelengths
fName = fullfile(isetbioDataPath, 'human', 'luminosity.mat');
V = ieReadSpectra(fName, wave);

irradianceP = oiGet(oi, 'photons');
if isempty(irradianceP)
    illuminance = [];
    meanIlluminance = [];
    return;
end

try
    % Formula requires irradiance in energy units
    irradianceE = Quanta2Energy(wave, irradianceP);

    % Do the calculation.
    img = RGB2XWFormat(irradianceE);
    illuminance = (683 * binWidth) * img * V;
    illuminance = XW2RGBFormat(illuminance, sz(1), sz(2));
catch
    % We are probably here because of a memory problem. So, let's try the
    % calculation again, but one waveband at a time
    [r, c, w] = size(irradianceP);
    illuminance = zeros(r, c);
    clear irradianceP;

    for ii = 1:w
        irradianceP = oiGet(oi, 'photons', wave(ii));
        illuminance = illuminance + (683 * binWidth) * ...
            Quanta2Energy(wave(ii), irradianceP) * V(ii);
    end
end

% Compute the mean if requested.
if nargout > 1, meanIlluminance = mean(illuminance(:)); end

% Compute the complementary (infrared mainly) illuminance if requested
if nargout > 2
    % shiftedV = V;
    % oldPeak = find(V == max(V)); % The luminosity function's peak
    newPeak = find(wave > 750); % Move peak to the right to here
    if isempty(newPeak)
        return;
    else
        % rightShift = newPeak(1) - oldPeak;
        % shiftedV = circshift(V, rightShift);

        try
            % Formula requires irradiance in energy units
            irradianceE = Quanta2Energy(wave, irradianceP);

            % Do the calculation.
            img = RGB2XWFormat(irradianceE);

            % Invisible energy
            % compIlluminance = (683 * binWidth) * img * (1 - V);

            % Shifted luminosity
            % compIlluminance = (683 * binWidth) * img * (shiftedV);

            % Total energy
            compIlluminance = (683 * binWidth) * img;

            compIlluminance = XW2RGBFormat(compIlluminance, sz(1), sz(2));

        catch
            % We are probably here because of a memory problem.  We should
            % check the Matlab Error. At this point, we simply assume so
            % and then we try the calculation again, but one waveband at a
            % time
            [r, c, w] = size(irradianceP);
            compIlluminance = zeros(r, c);
            clear irradianceP;

            %% m.p. 12/16/2007
            % We need to plug in a function here that describes a new
            % metric for quantifying the effect of IR energy.

            for ii = 1:w
                irradianceP = oiGet(oi, 'photons', wave(ii));

                % Invisible energy
                % compIlluminance = compIlluminance + (683 * binWidth) ...
                %    * Quanta2Energy(wave(ii), irradianceP .* (1 - V(ii)));

                % Total energy
                compIlluminance = compIlluminance + (683 * binWidth) * ...
                    Quanta2Energy(wave(ii), irradianceP * (1 - V(ii)));

                % Shifted luminosity
                % compIlluminance = compIlluminance + (683 * binWidth) ...
                %    * Quanta2Energy(wave(ii), irradianceP) * ...
                %    (1 - shiftedV(ii));
            end
        end
    end
    meanCompIlluminance = mean(compIlluminance(:));
end

end
