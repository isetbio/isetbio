function spdPerWlBand = rtbSpdPowerPerNmToPowerPerWlBand(spdPerNm, S)
%% Convert spectral power-per-nanometer to power-per-wavelength-band
%
% spdPerWlBand = rtbSpdPowerPerNmToPowerPerWlBand(spdPerNm, S)
% Converts the given spdPerNm matrix, which should contain a spectral
% power distribution with samples in units of power-per-nanometer, to the
% equivalent distribution with in units of power-per-wavelength-band.  The
% given S must describe the spectral sampling used in spdPerNm, and
% determines the correct conversion factor.
%
% Returns the given spdPerNm, multiplied with the spectral sampling band
% width in the given S.
%
% spdPerWlBand = rtbSpdPowerPerNmToPowerPerWlBand(spdPerNm, S)
%
%%% RenderToolbox3 Copyright (c) 2012-2013 The RenderToolbox3 Team.
%%% About Us://github.com/DavidBrainard/RenderToolbox3/wiki/About-Us
%%% RenderToolbox3 is released under the MIT License.  See LICENSE.txt.

parser = inputParser();
parser.addRequired('spdPerNm', @isnumeric);
parser.addRequired('S', @isnumeric);
parser.parse(spdPerNm, S);
spdPerNm = parser.Results.spdPerNm;
S = parser.Results.S;

% get the sampling bandwidth
S = MakeItS(S);
bandwidth = S(2);

% multiply band power with bandwidth to get power per spectrum band
spdPerWlBand = spdPerNm  .* bandwidth;
