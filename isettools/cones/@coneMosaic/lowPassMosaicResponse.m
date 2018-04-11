function [lowPassedResponse, Lmap, Mmap, Smap] = ...
    lowPassMosaicResponse(obj, absorptions, spaceConstants)
% Low pass filter mosaic isomerizations
%
% Syntax:
%    [lowPassedResponse, Lmap, Mmap, Smap] = ...
%         lowPassMosaicResponse(obj, response, spaceConstants)
%
% Description:
%    Spatially low pass filter a mosaic's response using different space
%    constants for each of the L, M, and S-cone submosaic.  This is done
%    after demosaicing the absorptions using the
%    coneMosaic.demosaiedResponses method.
%
%    A re-mosaiced version of the low-passed maps is returend, as well as
%    low-passed demosaiced images corresponding to the L, M and S cones.
%
% Inputs:
%    obj               - A coneMosaic object
%    absorptions       - The absorptions matrix from the mosaic.
%                        Absorptions is a synonym for isomerizations in
%                        ISETBio.
%    spaceConstants    - Three vector of space constants for L, M, and S
%                        cones in microns.
%
% Outputs:
%    lowPassedResponse - Remosaiced low-passed absorptions 
%    Lmap              - Demosaiced low-passed L cone plane  
%    Mmap              - Demosaiced low-passed M cone plane  
%    Smap              - Demosaiced low-passed S cone plane    
%
% Optional key/value pairs:
%     None.
%
% Notes:
%    * [NOTE: DHB - Need to say what kind of filter is applied and what the
%      space constant parameter means in the context of the filter
%      parameterization.]
%    * [NOTE: DHB - Please say what kind of function this is.] (from below)

% History:
%    xx/xx/16  NPC  ISETBIO Team, 2016
%    08/08/17  dhb  Figured out what the inputs & outputs are and commented
%    02/23/18  jnm  Formatting

    %% Get desmosaiced responses
    demosaicedResponses = obj.demosaicedResponses(absorptions);

    %% Initialize
    lowPassedResponse = 0 * absorptions;
    spaceAxis = obj.patternSupport(1, :, 1) * 1e6;
    responseSize = size(absorptions);

    %% Loop over cone types
    for mosaicIndex = 1:3
        % Create low pass filter
        lowPassFilter = generateLowPassFilter(...
            ceil(spaceConstants(mosaicIndex)) / ...
            (spaceAxis(2) - spaceAxis(1)));

        % Pad the demosaiced responses for this cone type
        margin = size(lowPassFilter, 2);
        paddedResponse = padarray(squeeze(...
            demosaicedResponses(:, :, mosaicIndex)), margin*[1 1], 0);

        % Convolve
        tmp = conv2(paddedResponse, lowPassFilter, 'same');
        tmp = tmp(margin + (1:responseSize(1)), ...
            margin + (1:responseSize(2)));

        % Pop the result into the places where there is a cone of this type
        indices = find(obj.pattern == mosaicIndex + 1);
        lowPassedResponse(indices) = tmp(indices);

        % Set returned demosaiced low-pass filtered maps for each cone type
        if (mosaicIndex == 1)
            Lmap = tmp;
            Lmap(1:size(lowPassFilter, 1), 1:size(lowPassFilter, 2)) = ...
                lowPassFilter / max(lowPassFilter(:)) * max(Lmap(:));
        elseif (mosaicIndex == 2)
            Mmap = tmp;
            Mmap(1:size(lowPassFilter, 1), 1:size(lowPassFilter, 2)) = ...
                lowPassFilter / max(lowPassFilter(:)) * max(Mmap(:));
        else 
            Smap = tmp;
            Smap(1:size(lowPassFilter, 1), 1:size(lowPassFilter, 2)) = ...
                lowPassFilter / max(lowPassFilter(:)) * max(Smap(:));
        end
    end

end

%% Generate the convolution kernel
%
% [DHB NOTE: Please say what kind of function this is.]
function kernel2D = generateLowPassFilter(sigma)
% generate the low pass filter
%
% Syntax:
%   kernel2D = generateLowPassFileter(sigma)
%
% Description:
%    Generate the convolution kernel by generating a low pass filter
%
% Inputs:
%    sigma    - The modifier
%
% Outputs:
%    kernel2D - The generated 2D kernel
%
% Optional key/value pairs:
%    None.
%
    s = sqrt(-log(0.001) / 0.5);   % 1/1000
    spaceAxis = -round(s * sigma):1:round(s * sigma);
    kernel = exp(-0.5 * (spaceAxis / sigma) .^ 2);
    kernel2D = kernel' * kernel;
    kernel2D = kernel2D / sum(kernel2D(:));
end