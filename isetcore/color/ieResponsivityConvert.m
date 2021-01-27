function [responsivity, sFactor] = ...
    ieResponsivityConvert(responsivity, wave, method)
% Convert sensory responsivity in photons to energy, or vice versa
%
% Syntax:
%   [responsivity, sFactor] = ...
%       ieResponsivityConvert(responsivity, wave, [method]);
%
% Description:
%    When calculating a sensor responsivity, it is essential to specify
%    whether the input is in units of photons or energy. In a digital
%    imager, for example, each photon has a probability of producing an
%    electron. Different wavelengths must account for the number of
%    photons, even if the data are expressed in units of energy.
%
%    Isetbio uses photons as the basis for nearly all of the response
%    calculations.
%
%    But some important sensors are defined with respect to signal energy.
%    The most important of these are the XYZ sensors. These are specified
%    with respect to energy. It is also the case that the human cone
%    responses are often specified with respect to energy units.
%
%    In some cases in the code, we convert the input signal in photons to
%    energy and use the standard XYZ values.
%
%    In other cases, however, we have many inputs and it is easier to
%    convert the specification of the XYZ functions into a form that is
%    correct for photon inputs. This function performs that conversion.
%
%    Suppose transQ, transE are responsivity measured with respect to
%    quanta and energy. Suppose that E2Q is the conversion from energy to
%    quanta as a function of wavelength. Finally, suppose inE and inQ are
%    input signals in energy and quanta units.
%
%     response = transE' * inE = (transE' * (1 / E2Q)) * (E2Q * inE) ...
%              = transQ' * inQ
%
%    We can see that transQ is related to transE as
%       transQ' = transE' * (1 / E2Q).
%
%    This routine converts responsivities measured in energy units (respE)
%    to responsivities appropriate for photons calculations (respQ).
%
%    These issues are handled explicitly in ieLuminanceFromEnergy,
%    ieLuminanceFromPhotons and ieXYZFromEnergy
%
%    To specify filter transmissivities, it is not necessary to pay
%    attention to the input signal units (photons or energy). Filters
%    transmit a fraction of the photons and they transmit the same fraction
%    of the energy.
%
%    The color responsivities are in the columns of the variable
%    responsivity. The variable wave is the wavelength in nanometers. If
%    the variable method = 'e2q' this routine converts filters specified
%    for energy to work with photons (quanta). if method = 'q2e' this
%    routine converts filters for quanta to work with energy.
%
%   This function contains examples of usage inline. To access these, enter
%   'edit ieResponsivityConvert.m' into the Command Window.
%
% Inputs:
%    responsivity - Matrix. The columns represent color responsivities.
%    wave         - Vector. A wavelength sampling (nm) as a vector.
%    method       - (Optional) String. A representation of the desired
%                   methods. The options include:
%        'e2q': (Default) Filters specified for energy to work with
%                quanta/photons
%        'q2e': Filters specified for quanta to work with energy
%
% Outputs:
%    responsivity - Matrix. The columns of color responsivities after
%                   manipulation to the corresponding new format
%    sFactor      - Vector. The sensitivity factor. A row vector of how
%                   how the input was scaled into the output at each
%                   wavelength, before a final overall scaling.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note - DHB: There is a rescaling of the output at the end of this
%      routine that makes the max value of the output have the same value
%      as the max value of the input sensitivities. This makes it so that
%      computing with the input and output sensitivies does not produce
%      the same answer. I don't think this scaling is a good idea, but
%      left it alone in case other code counts on it. Currently, the
%      example at the top of the source code, which is meant to show that
%      you get the same answer working in either energy or quantal units
%      when you do the conversions right, does not show this because of
%      the rescaling.]
%    * [Note - DHB: To make the example work, I had to read in some energy
%      sensitivities, the old example method used a function that no
%      longer exists. I think ieReadSpectral('stockman',wave) gets cone
%      sensitivities in energy units, but there doesn't appear to be any
%      easy way to find out. Check and fix if necessary.]
%
% See Also:
%   ieLuminanceFromEnergy, ieLuminanceFromPhotons, ieXYZFromEnergy.
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC.
%    10/27/17  jnm  Comments & formatting
%    11/16/17  jnm  Formatting
%    07/11/19  JNM  Formatting update

% Examples:
%{
   % Signal in photon units
   wave = 400:10:700;
   signalEnergy = ieReadSpectra('D65', wave);
   signalPhotons = Energy2Quanta(wave(:), signalEnergy(:));

   % Signal in energy units
   conesE = ieReadSpectra('stockman', wave);
   [conesP, sFactor] = ieResponsivityConvert(conesE, wave, 'e2q');

   % These two calculations produce equal results
   vP = conesP' * signalPhotons(:)
   vE = conesE' * signalEnergy(:)
%}

if notDefined('responsivity')
    error('Must define color responsivity functions');
end
if notDefined('wave')
    error('Must define wavelength sampling in nanometers');
end
if notDefined('method'), method = 'e2q'; end

if length(wave) ~= size(responsivity, 1)
    error('Mis-match between wavelength and color filters.');
end

maxTrans = max(responsivity(:));
switch lower(method)
    case {'e2q', 'energy2quanta', 'e2p', 'energy2photons'}
        % Set up filters that handle energy to handle quanta
        sFactor = Quanta2Energy(wave(:), ones(1, length(wave)));
    case {'q2e', 'quanta2energy', 'p2e', 'photons2energy'}
        % Set up filters that handle energy to handle quanta
        sFactor = Energy2Quanta(wave(:), ones(1, length(wave))');
    otherwise
        error('Unknown method');
end
responsivity = diag(sFactor) * responsivity;

% The throughput at max should be the same
responsivity = responsivity * (maxTrans / (max(responsivity(:))));

end