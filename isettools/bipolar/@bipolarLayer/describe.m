function str = describe(obj)
% Summarize the bipolar layer properties in a string% Syntax:
%
% Syntax:
%   @bipolarLayer.describe;
%
% Description:
%    Summarize the bipolar layer properties in a string
%
% Inputs:
%   None
%
% Outputs:
%    str - A string containing the patch size, center, time step, number of
%          trials, eye side, species, and mosaic-specific row & column
%          properties of the object
% Notes:
% * Do we want to have the option to include more parameters (or fewer) as
%   options for describe?

% History: 
% BW ISETBIO Team, 2017
%
%    10/18/17  jnm  Comments & Formatting

%%
gdata = guidata(obj.fig);
nMosaic = get(gdata.listMosaics, 'Value');

%% General properties
str = [];
txt = sprintf('Patch size:   \t%.1f %.1f um\n', obj.input.size * 1e6);
str = addText(str, txt);
txt = sprintf('Center:       \t(%.1f, %.1f) mm\n', ...
    [obj.input.center] * 1e3);
str = addText(str, txt);
txt = sprintf('Time step:    \t%.1f ms\n', ...
    obj.input.integrationTime * 1e3);
str = addText(str, txt);
txt = sprintf('N Trials:     \t%.0f \n', obj.nTrials);
str = addText(str, txt);
txt = sprintf('Eye side:     \t%s \n', obj.input.whichEye);
str = addText(str, txt);
txt = sprintf('Species:      \t%s \n', obj.input.species);
str = addText(str, txt);

%% Mosaic-specific properties
txt = sprintf('---\n');
str = addText(str, txt);
[r, c, ~] = size(obj.mosaic{nMosaic}.responseCenter);
txt = sprintf('Row, Col:   \t%d %d \n', r, c);
str = addText(str, txt);

end
