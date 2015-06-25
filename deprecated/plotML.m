function plotML(ml,type)
%
%   plotML(ml,type)
%
% Author: ImagEval
% Purpose:
%   Various microlens summary plots
%
% Example
%   plotML(ml,'pixelIrradiance');


% Make a figure showing the etendue across the array.  The units of the
% support are unclear to me at this moment.
figNum = vcSelectFigure('GRAPHWIN');
plotSetUpWindow(figNum);

switch lower(type)
    case {'pixelirradiance'}
        irrad = mlensGet(ml,'pixelIrradiance');
        x = mlensGet(ml,'xcoordinate');
        mesh(x,x,irrad);
        xlabel('Position (um)');
        ylabel('Position (um)');
        zlabel('Relative irradiance');
        h = hot(256); colormap(h(30:220,:))
    otherwise
        error('Unknown plotML type.');
end

return;
