% Method that low-passes a mosaic's response using different space
% constants for each of the L,M, and S-cone submosaic.
% 
% NPC, ISETBIO Team, 2016  
%

function [lowPassedResponse,  Lmap, Mmap, Smap] = lowPassMosaicResponse(obj, response, spaceConstants)
    
    demosaicedResponses = obj.demosaicedResponses(response);
    lowPassedResponse = 0*response;
    spaceAxis = obj.patternSupport(1,:,1)*1e6;
    responseSize = size(response);
    
    for mosaicIndex = 1:3
        lowPassFilter = generateLowPassFilter(ceil(spaceConstants(mosaicIndex)) / (spaceAxis(2)-spaceAxis(1)));
        margin = size(lowPassFilter,2);
        paddedResponse = padarray(squeeze(demosaicedResponses(:,:,mosaicIndex)), margin*[1 1], 0);
        tmp = conv2(paddedResponse, lowPassFilter, 'same');
        tmp = tmp(margin+(1:responseSize(1)), margin+(1:responseSize(2)));
        indices = find(obj.pattern == mosaicIndex+1);
        lowPassedResponse(indices) = tmp(indices);
        if (mosaicIndex == 1)
            Lmap = tmp;
            Lmap(1:size(lowPassFilter,1), 1:size(lowPassFilter,2)) = lowPassFilter/max(lowPassFilter(:)) * max(Lmap(:));
        elseif (mosaicIndex == 2)
            Mmap = tmp;
            Mmap(1:size(lowPassFilter,1), 1:size(lowPassFilter,2)) = lowPassFilter/max(lowPassFilter(:)) * max(Mmap(:));
        else 
            Smap = tmp;
            Smap(1:size(lowPassFilter,1), 1:size(lowPassFilter,2)) = lowPassFilter/max(lowPassFilter(:)) * max(Smap(:));
        end
    end
    
end

function kernel2D = generateLowPassFilter(sigma)
    s = sqrt(-log(0.001)/0.5);   % 1/1000
    spaceAxis   = -round(s*sigma):1:round(s*sigma);
    kernel      = exp(-0.5*(spaceAxis/sigma).^2);
    kernel2D    = kernel' * kernel;
    kernel2D    = kernel2D / sum(kernel2D(:));
end