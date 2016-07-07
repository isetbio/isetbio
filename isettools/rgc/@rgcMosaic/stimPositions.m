function [stimX, stimY, offset] = stimPositions(rgcMosaic,xcell,ycell)
% Calculate the receptive field positions in the input frame for one RGC
%
%  [stimX, stimY, offset] = stimPositions(obj,xcell,ycell)
%
% @JRG More comments please.
%
% JRG ISETBIO Team

% The RGC center location
stimCenterCoords = rgcMosaic.cellLocation{xcell,ycell};

% Find the spatial extent of the RF in terms of multiples of rfDiameter
extent = 1; % Set spatial RF extent

% Pull out stimulus coordinates of interest

% Get midpoint of RF by taking half of the col size
sRFMidPointX = floor((extent/2)*size(rgcMosaic.sRFcenter{1,1},1));

% stimCenterCoords indicates position of RGC on stimulus image
% The first x coord of the stimulus of interest is the RGC center minus
% the midpoint size of the RF.
xStartCoord = (stimCenterCoords(1) - sRFMidPointX);
xEndCoord   = (stimCenterCoords(1) + sRFMidPointX);

% Get midpoint of RF by taking half of the row size
sRFMidPointY = floor((extent/2)*size(rgcMosaic.sRFcenter{1,1},2));

% stimCenterCoords indicates position of RGC on stimulus image
% The first x coord of the stimulus of interest is the RGC center minus
% the midpoint size of the RF.
yStartCoord = (stimCenterCoords(2) - sRFMidPointY);
yEndCoord   = (stimCenterCoords(2) + sRFMidPointY);

stimX =  floor(xStartCoord):floor(xEndCoord);
stimY =  floor(yStartCoord):floor(yEndCoord);

if nargout == 3
    % Explain the offset better
    % @JRG
    
    % Set rounding of cell location based on whether it is
    % positive or negative
    if rgcMosaic.cellLocation{1,1}(1) > 0
        offset(1) = ceil(rgcMosaic.cellLocation{1,1}(1));
    else
        offset(1) = floor(rgcMosaic.cellLocation{1,1}(1));
    end
    
    if rgcMosaic.cellLocation{1,1}(2) > 0
        offset(2) = ceil(rgcMosaic.cellLocation{1,1}(2));
    else
        offset(2) = floor(rgcMosaic.cellLocation{1,1}(2));
    end
end

end