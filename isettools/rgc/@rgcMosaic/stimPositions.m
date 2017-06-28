function [stimX, stimY, offset] = stimPositions(rgcMosaic,xcell,ycell)
% Calculate the receptive field positions in the input frame for one RGC
%
%  [stimX, stimY, offset] = stimPositions(obj,xcell,ycell)
%
% To compute the linear response of an RGC, the relevant spatial positions
% of the stimulus must be determined based on the center coordinates of the
% RF in stimulus space and the spatial extent of the RF. This function
% pulls out the relevant spatial coordinates for a given RGC.
%
% 5/2016 JRG (c) ISETBIO Team

% If we want the stimulus center in units of microns, I am worried.
% micronsToBipolars has units of cell/micron.  cell Location has units of
% cell.  We are multiplying them.  (BW)
%
% Also, the logic below around if/else and ceil() seems broken to me. (BW).

% This code is very confusing to me (BW) and I hope to straighten out the
% logic with JRG.
% The RGC center location
micronsToBipolars = size(rgcMosaic.cellLocation,2)/(1e6*rgcMosaic.Parent.size(2));
stimCenterCoords = micronsToBipolars*rgcMosaic.cellLocation{xcell,ycell};

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

stimX =  ceil(xStartCoord):floor(xEndCoord);
stimY =  ceil(yStartCoord):floor(yEndCoord);

if length(stimX)>length(stimY); stimX = stimX(1:length(stimY)); end
if length(stimY)>length(stimX); stimY = stimY(1:length(stimX)); end

if length(stimX)>size(rgcMosaic.sRFcenter{xcell,ycell},1) || length(stimY)>size(rgcMosaic.sRFcenter{xcell,ycell},2) 
    stimX = stimX(1:size(rgcMosaic.sRFcenter{xcell,ycell},1)); 
    stimY = stimY(1:size(rgcMosaic.sRFcenter{xcell,ycell},2)); 
end
% if length(stimY)>size(rgcMosaic.sRFcenter{xcell,ycell},2); stimY = stimY(1:size(rgcMosaic.sRFcenter{xcell,ycell},2)); end;

if nargout == 3
    % An offset is sometimes used because RGC mosaics may be defined with
    % their center coordinates not at (0,0).
    
    % Set rounding of cell location based on whether it is
    % positive or negative
    % BW.  Hunh?  Aren't the functions the same on the if/else?

    if rgcMosaic.cellLocation{1,1}(1) > 0
        offset(1) = ceil(rgcMosaic.cellLocation{1,1}(1));
    else
        offset(1) = ceil(rgcMosaic.cellLocation{1,1}(1));
    end
    
    if rgcMosaic.cellLocation{1,1}(2) > 0
        offset(2) = ceil(rgcMosaic.cellLocation{1,1}(2));
    else
        offset(2) = ceil(rgcMosaic.cellLocation{1,1}(2));
    end
end

end