function [wRange nRange cutoff] = rgcFindGoodDistanceParameters(layer,assign)
% Tries to find good parameters for the distance function
%
%   [wRange nRange cutoff] = rgcFindGoodDistanceParameters(layer,[assign])
%
% The distance function XXX comments here ...
%
% layer:  rgcLayerObject
% assign: do we assign the values after computation?
%         default = 0 = 'no'
%
% If the function raises 'Cannot find good parameters' error, try new
% starting parameters
%
% (c) 2010 Stanford Synapse Team

if (notDefined('layer') || ~isequal(class(layer),'rgcLayer'))
    error('rgcLayer class required.');
end
if notDefined('assign'),  assign = 0; end
if isequal('yes',assign), assign = 1; end
if isequal('no',assign),  assign = 0; end

wRange  = layer.get('wRange');
nRange  = layer.get('nRange');
cutoff = layer.get('cutoff'); 
distanceFunction = layer.parent.get('distanceFunction');
cellSpacing = layer.get('cellSpacing');

cubeSide=0;
threshold = cutoff-nRange;
while distanceFunction(cubeSide*cellSpacing,wRange)>=threshold
    cubeSide = cubeSide+1;
end
cubeSide = max(cubeSide-1,0);
    
previousSide=cubeSide;
previousc = cutoff;
previousw = wRange;
previousn = nRange;

nTest = 0;

while cubeSide ~= 2
    nTest = nTest + 1;
    if nTest == 1000
        error('Cannot find good parameters');
    end
    
    if (previousSide > 2 && cubeSide > 2)
        cutoff = 2*cutoff;
        nRange = nRange/2;
        wRange = wRange/2;
    elseif (previousSide < 2 && cubeSide < 2)
        cutoff = cutoff/2;
        nRange = nRange*2;
        wRange = wRange*2;
    else
        cutoff = (cutoff+previousc)/2;
        nRange = (nRange+previousn)/2;
        wRange = (wRange+previousw)/2;
    end
    if cutoff/2 < nRange
        nRange = cutoff/2;
    end
    
    previousSide = cubeSide;
    cubeSide=0;
    threshold = cutoff-nRange;
    
    while distanceFunction(cubeSide*cellSpacing,wRange)>=threshold
        cubeSide = cubeSide+1;
    end
    cubeSide = max(cubeSide-1,0);
end

if assign
    layer.set('wRange',wRange);
    layer.set('nRange',nRange);
    layer.set('cutoff',cutoff);
end

return