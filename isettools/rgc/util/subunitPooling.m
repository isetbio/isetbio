function [spResponseCenter, spResponseSurround] = subunitPooling(spResponseCenterLinear, spResponseSurroundLinear, varargin)
% Nonlinear spatial summation of subunit receptive fields. The user may
% choose the subunit model to implement.
%
% Subunit options:
% 1. 'pixel': Every pixel in spatial response goes through nonlinearity before
% summation, as in Meister & Gollisch 2008, and Baccus, Olveiczky and
% Meister 2008.
%
% 2. 'surround': linear spatial center, sRFsurround is pooled into two subunits on the
% left and right, Takeshita & Gollisch 2014.
%
% 3/2016 JRG (c) isetbio team

p = inputParser;
p.addRequired('spResponseCenterLinear');
p.addRequired('spResponseCenterSurround');
p.addParameter('model', 'pixel', @ischar);
p.parse(spResponseCenterLinear, spResponseSurroundLinear, varargin{:});
model = p.Results.model;

switch ieParamFormat(model)
    case{'pixel'}
        % Every pixel in spatial response goes through nonlinearity before summation,
        % Meister & Gollisch 2008, and Baccus, Olveiczky and Meister 2008

        spResponseCenter   = cellfun(@(x) exp(.014*x), spResponseCenterLinear,'uniformoutput',false);
        spResponseSurround = cellfun(@(x) exp(.014*x), spResponseSurroundLinear,'uniformoutput',false);
        
    case{'surround'}
        % linear spatial center, sRFsurround is pooled into two subunits on the left and right
        % Takeshita & Gollisch 2014        
        
        nCells = size(spResponseCenterLinear);
        spResponseSize = size(spResponseSurroundLinear{1,1});
        for xc = 1:nCells(1)
            for yc = 1:nCells(2)
                
                spResponseCenterSub{xc,yc} = squeeze(exp(.014*mean(mean(spResponseCenterLinear{xc,yc},1),2)));
                
                spResponseSurroundSub{xc,yc} = squeeze(...
                    sum(sum(exp(.014*spResponseSurroundLinear{xc,yc}(:,1:floor(spResponseSize(2)/2),:,:)),1),2) + ...
                    sum(sum(exp(.014*spResponseSurroundLinear{xc,yc}(:,ceil(spResponseSize(2)/2):spResponseSize(2),:,:)),1),2));
                %                         spResponseFull = spResponseCenter{xc,yc} + spResponseSurround{xc,yc};
                %                         spResponseVec{xc,yc} = squeeze(mean(mean(spResponseFull,1),2))';
                
                for frameNumber = 1:spResponseSize(3)
                    for rgbInd = 1:3
                        spResponseCenter{xc,yc}(:,:,frameNumber,rgbInd) = ...
                            (spResponseCenterSub{xc,yc}(frameNumber,rgbInd)/(spResponseSize(1)*spResponseSize(2)))*...
                            ones([spResponseSize(1) spResponseSize(2) 1 1]);
                        spResponseSurround{xc,yc}(:,:,frameNumber,rgbInd) = ...
                            (spResponseSurroundSub{xc,yc}(frameNumber,rgbInd)/(spResponseSize(1)*spResponseSize(2)))*...
                            ones([spResponseSize(1) spResponseSize(2) 1 1]);
                        % permute(repmat(spResponseSurroundSub{xc,yc},[1 1 spResponseSize(1) spResponseSize(2)]),[3 4 1 2]);
                    end
                end
                
            end
        end
        % Need to fix center response
        spResponseCenter   = cellfun(@(x) exp(.014*x), spResponseCenterLinear,'uniformoutput',false);
                
    case{'surroundfour'}
        % After Shah & Chichilnisky, manuscript      
        
        nSubunitsX = 2; nSubunitsY = 2;
        
        nCells = size(spResponseCenterLinear);
        spResponseSize = size(spResponseSurroundLinear{1,1});
        
        pixelsPerSubunitRow = spResponseSize(1)/nSubunitsX;
        pixelsPerSubunitCol = spResponseSize(2)/nSubunitsY;
        
        for xc = 1:nCells(1)
            for yc = 1:nCells(2)
                
                spResponseCenterSub{xc,yc} = squeeze(exp(mean(mean(spResponseCenterLinear{xc,yc},1),2)));
                
                spResponseSurroundSub{xc,yc} = 0;
                for iSubunitX = 1:nSubunitsX
                    for iSubunitY = 1:nSubunitsY
                        
                        xcoords = ceil((iSubunitX-1)*pixelsPerSubunitRow + 1) : floor((iSubunitX)*pixelsPerSubunitRow);
                        ycoords = ceil((iSubunitY-1)*pixelsPerSubunitCol + 1) : floor((iSubunitY)*pixelsPerSubunitCol);
                        
                        spResponseSurroundSub{xc,yc} = spResponseSurroundSub{xc,yc} + squeeze(...
                            sum(sum(exp(.014*spResponseSurroundLinear{xc,yc}(xcoords,ycoords,:,:)),1),2));
                    end
                end
                
                for frameNumber = 1:spResponseSize(3)
                    for rgbInd = 1:3
                        spResponseCenter{xc,yc}(:,:,frameNumber,rgbInd) = (...
                            (spResponseCenterSub{xc,yc}(frameNumber,rgbInd)/(spResponseSize(1)*spResponseSize(2)))*...
                            ones([spResponseSize(1) spResponseSize(2) 1 1]));
                        spResponseSurround{xc,yc}(:,:,frameNumber,rgbInd) = (...
                            (spResponseSurroundSub{xc,yc}(frameNumber,rgbInd)/(spResponseSize(1)*spResponseSize(2)))*...
                            ones([spResponseSize(1) spResponseSize(2) 1 1]));
                        
                    end
                end
                
            end
        end
        % Need to fix center response
        spResponseCenter   = cellfun(@(x) exp(.014*x), spResponseCenterLinear,'uniformoutput',false);
end