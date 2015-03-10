function plotMLOffsets(ISA)
%
%   plotMLOffsets(ISA)
%
% Plot a mesh graph showing the optimal offsets of the microlens
% positions.  Negative means displaced towards the center.
%
%
% Copyright ImagEval Consultants, LLC, 2005.

if notDefined('ISA'), ISA = vcGetObject('ISA'); end

ml = sensorGet(ISA,'ml');
optimalOffsets = mlensGet(ml,'microoptimaloffsets');

figNum = vcSelectFigure('GRAPHWIN');
plotSetUpWindow(figNum);

support = sensorGet(ISA,'spatialSupport','mm');
mesh(support.y, support.x, optimalOffsets);
xlabel('Position (mm)');
ylabel('Position (mm)');
zlabel('Optimal Microlens Offset (um)');

%
uData.support = support;
uData.optimalOffsets = optimalOffsets;
uData.command = 'mesh(support.y, support.x, optimalOffsets)';

set(figNum,'userdata',uData);

return;
