%%% RenderToolbox3 Copyright (c) 2012-2013 The RenderToolbox3 Team.
%%% About Us://github.com/DavidBrainard/RenderToolbox3/wiki/About-Us
%%% RenderToolbox3 is released under the MIT License.  See LICENSE.txt.
%
% Convert multi-spectral data to a sensor image representation.
%   @param multispectralData multispectral image matrix
%   @param imageS image spectral sampling description
%   @param matchingFunction color matching funciton matrix or filename
%   @param matchingS matching function spectral sampling description
%
% @details
% Convert the given @a multispectralData to a sensor image representation,
% based on the given @a matchingFunction.  @a multispectralData should
% have size [height width nImageSpecralSamples] and spectral sampling
% described by @a imageS of the form [start delta nImageSpecralSamples].
% @a multispectralData should contain data in units of Power per Unit
% Wavelength.
%
% @details
% @a matchingFunction may have one of two forms.  It may be a matrix with
% size [nChannels, nMatchingSpecralSamples].  When @a matchingFunction
% is a matrix, @a matchingS must be provided as a description of the
% matching function spectral sampling, in the form [start delta
% nMatchingSpecralSamples].
%
% @details
% @a matchingFunction may also be the name of a Psychtoolbox colorimetric
% data file that contains a matrix of matching function data and an "S"
% descripton of the matching function spectral sampling.  By convention,
% these data files are named with the prefix "T_".  When @a
% matchingFunction is a file name, @a matchingS is ignored.
%
% @details
% For more about Psychtooblox colorimetric .mat files and conventions, see
% the <a
% href="http://docs.psychtoolbox.org/PsychColorimetricMatFiles">Psychtoolbox web documentation</a>
% or the file
%   Psychtoolbox/PsychColorimetricData/PsychColorimetricMatFiles/Contents.m
%
% @details
% Returns a sensor image representation of the given @a multispectralData,
% with size [height width nSensorChannels].
%
% @details
% Usage:
%   sensorImage = MultispectralToSensorImage(multispectralData, imageS, matchingFunction, matchingS)
%
% @ingroup Utilities
function sensorImage = MultispectralToSensorImage(multispectralData, imageS, matchingFunction, matchingS)
% load matchingFunction and matchingS from Psychtoolbox .mat file?
if ischar(matchingFunction)
    [matchingFunction, matchingS] = ...
        ParsePsychColorimetricMatFile(matchingFunction);
end

% "multiply in" to Psychtoolbox convention of Power per Wavelength Band
multispectralData = ...
    SpdPowerPerNmToPowerPerWlBand(multispectralData, imageS);

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
