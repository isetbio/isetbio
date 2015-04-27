%%% RenderToolbox3 Copyright (c) 2012-2013 The RenderToolbox3 Team.
%%% About Us://github.com/DavidBrainard/RenderToolbox3/wiki/About-Us
%%% RenderToolbox3 is released under the MIT License.  See LICENSE.txt.
%
% Convert spectral power-per-nanometer to power-per-wavelength-band
%   @param spdPerNm matrix of data with power-per-nanometer
%   @param S description of wavelength sampling, [start delta n]
%
% @details
% Converts the given @a spdPerNm matrix, which should contain a spectral
% power distribution with samples in units of power-per-nanometer, to the
% equivalent distribution with in units of power-per-wavelength-band.  The
% given @a S must describe the spectral sampling used in @a spdPerNm, and
% determines the correct conversion factor.
%
% @details
% Returns the given @a spdPerNm, multiplied with the spectral sampling band
% width in the given @a S.
%
% @details
% Usage:
%   spdPerWlBand = SpdPowerPerNmToPowerPerWlBand(spdPerNm, S)
%
% @ingroup Utilities
function spdPerWlBand = SpdPowerPerNmToPowerPerWlBand(spdPerNm, S)

% get the sampling bandwidth
S = MakeItS(S);
bandwidth = S(2);

% multiply band power with bandwidth to get power per spectrum band
spdPerWlBand = spdPerNm  .* bandwidth;
