function [stimX, stimY, offset] = stimPositions(rgcMosaic,xcell,ycell)
% Calculate the receptive field positions in the input frame for one RGC
%
%  [stimX, stimY, offset] = stimPositions(obj,xcell,ycell)
%
% JRG ISETBIO Team

% The RGC center location
stimCenterCoords = rgcMosaic.cellLocation{xcell,ycell};

% Find the spatial extent of the RF in terms of multiples of rfDiameter
if isa(rgcMosaic, 'rgcPhys')
    % rgcPhys is the case when the receptive fields
    % (rgcMosaic) are fit from an experiment.
    
    % If all the RFs are the same, we can use {1,1}.
    % Otherwise ...
    extent = round(size(rgcMosaic.sRFcenter{1,1},1)/rgcMosaic.rfDiameter);
    
    % Pull out stimulus coordinates of interest
    % Words needed here to explain
    stimX =  ceil((stimCenterCoords(1) - floor((extent/2)*rgcMosaic.rfDiameter))):floor((stimCenterCoords(1) + floor((extent/2)*rgcMosaic.rfDiameter )));%
    stimY =  ceil((stimCenterCoords(2) - floor((extent/2)*rgcMosaic.rfDiameter))):floor((stimCenterCoords(2) + floor((extent/2)*rgcMosaic.rfDiameter )));%
    
else  % rgcGLM, rgcLinear
    
    extent = 1; % Set spatial RF extent
    
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
    
    % Pull out stimulus coordinates of interest
    stimX =  floor((stimCenterCoords(1) - floor((extent/2)*size(rgcMosaic.sRFcenter{1,1},1)))):floor((stimCenterCoords(1) + floor((extent/2)*size(rgcMosaic.sRFcenter{1,1},1) )));
    stimY =  floor((stimCenterCoords(2) - floor((extent/2)*size(rgcMosaic.sRFcenter{1,1},2)))):floor((stimCenterCoords(2) + floor((extent/2)*size(rgcMosaic.sRFcenter{1,1},2))));
    
end

end