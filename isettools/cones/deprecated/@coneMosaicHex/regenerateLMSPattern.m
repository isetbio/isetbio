function regenerateLMSPattern(obj, LMSdensity, varargin)
    % Regenerate the cone identities of the cone mosaic hex object 
    %
    % Syntax:
    %   regenerateLMSPattern(obj, LMSdensity, sConeMinDistanceFactor, sConeFreeRadiusMicrons)
    %
    % Description:
    %    Regenerate the LMS pattern for the provided cone mosaic hex using
    %    the passed LMSdensity and sConeMinDistanceFactor
    %
    % Inputs:
    %    obj                    - The cone mosaic hex object
    %    LMSdensity             - An 1x3 vector with desited L:M:S cone densities
    %    sConeMinDistanceFactor - Min distance factor for S-cones (usually
    %                             2.5, which results in 7% S-cones)
    %    sConeFreeRadiusMicrons - The foveola radius in microns which
    %                             contains no S-cones. Default: 45
    % Outputs:
    %    None.
    %
    % Optional key/value pairs:
    %    None.
    %
    
    p = inputParser;
    p.addParameter('sConeMinDistanceFactor', 2.5, @isnumeric);
    p.addParameter('sConeFreeRadiusMicrons', 45, @isnumeric);
    p.addParameter('zeroSconeDensity', false, @islogical);
    p.addParameter('visualizeRegeneratedMosaic', false, @islogical);
    p.parse(varargin{:});
    
    sConeMinDistanceFactor = p.Results.sConeMinDistanceFactor;
    sConeFreeRadiusMicrons = p.Results.sConeFreeRadiusMicrons;
    visualizeRegeneratedMosaic = p.Results.visualizeRegeneratedMosaic;
    
    LMSdensity = LMSdensity / sum(LMSdensity);
    
    idx = find(obj.pattern > 1);
    totalConesNum = numel(idx);
    newLconesNum = round(LMSdensity(1) * totalConesNum);
    newMconesNum = round(LMSdensity(2) * totalConesNum);
    newSconesNum = totalConesNum - newLconesNum - newMconesNum;
    
    obj.pattern = obj.pattern * 0;
    randomIndices = idx(randperm(totalConesNum));
    obj.pattern(randomIndices(1:newLconesNum)) = 2;
    obj.pattern(randomIndices(newLconesNum+(1:newMconesNum))) = 3;
    obj.pattern(randomIndices(newLconesNum+newMconesNum+(1:newSconesNum))) = 4;
    obj.reassignConeIdentities(...
                    'sConeMinDistanceFactor', sConeMinDistanceFactor, ...
                    'sConeFreeRadiusMicrons', sConeFreeRadiusMicrons);
    if (visualizeRegeneratedMosaic)
        obj.visualizeGrid();
    end
end

