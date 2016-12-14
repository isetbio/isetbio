function sensorImage = rtbMultispectralToSensorImage(multispectralData, imageS, matchingFunction, matchingS)
%% Convert multi-spectral data to a sensor image representation.
%
% sensorImage = rtbMultispectralToSensorImage(multispectralData, imageS, matchingFunction, matchingS)
% Converts the given multispectralData to a sensor image representation,
% based on the given matchingFunction.
%
% multispectralData must have size [height width nImageSpecralSamples] and
% spectral sampling described by imageS of the form [start delta
% nImageSpecralSamples].  multispectralData must contain data in units of
% Power per Unit Wavelength.
%
% matchingFunction may have one of two forms.  It may be a matrix with
% size [nChannels, nMatchingSpecralSamples].  When matchingFunction
% is a matrix, matchingS must be provided as a description of the
% matching function spectral sampling, in the form [start delta
% nMatchingSpecralSamples].
%
% matchingFunction may also be the name of a Psychtoolbox colorimetric
% data file that contains a matrix of matching function data and an "S"
% descripton of the matching function spectral sampling.  By convention,
% these data files are named with the prefix "T_".  When matchingFunction
% is a file name, matchingS is ignored.
%
% Returns a sensor image representation of the given @a multispectralData,
% with size [height width nSensorChannels].
%
% For more about Psychtooblox colorimetric .mat files and conventions, see
% the Psychtoolbox web documentation
%   http://docs.psychtoolbox.org/PsychColorimetricMatFiles
% or the file
%   Psychtoolbox/PsychColorimetricData/PsychColorimetricMatFiles/Contents.m
%
%%% RenderToolbox3 Copyright (c) 2012-2013 The RenderToolbox3 Team.
%%% About Us://github.com/DavidBrainard/RenderToolbox3/wiki/About-Us
%%% RenderToolbox3 is released under the MIT License.  See LICENSE.txt.

parser = inputParser();
parser.addRequired('multispectralData', @isnumeric);
parser.addRequired('imageS', @isnumeric);
parser.addRequired('matchingFunction');
parser.addRequired('matchingS', @isnumeric);
parser.parse(multispectralData, imageS, matchingFunction, matchingS);
multispectralData = parser.Results.multispectralData;
imageS = parser.Results.imageS;
matchingFunction = parser.Results.matchingFunction;
matchingS = parser.Results.matchingS;

% load matchingFunction and matchingS from Psychtoolbox .mat file?
if ischar(matchingFunction)
    [matchingFunction, matchingS] = ...
        rtbParsePsychColorimetricMatFile(matchingFunction);
end

% "multiply in" to Psychtoolbox convention of Power per Wavelength Band
multispectralData = ...
    rtbSpdPowerPerNmToPowerPerWlBand(multispectralData, imageS);

%% Convert image to arbitrary sensor representation by weighting.
% resample mathing function to match the image spectral sampling
imageS = MakeItS(imageS);
matchingS = MakeItS(matchingS);
matchingResampled = SplineCmf(matchingS, matchingFunction, imageS);

% weight multispectral planes by the matching function,
%   for each sensor channel
imageSize = size(multispectralData);
nChannels = size(matchingResampled, 1);
sensorImage = zeros(imageSize(1), imageSize(2), nChannels);
for ww = 1:imageS(3)
    for jj = 1:nChannels
        sensorImage(:,:,jj) = sensorImage(:,:,jj) ...
            + matchingResampled(jj,ww)*multispectralData(:,:,ww);
    end
end
