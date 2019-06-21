function str = describe(obj)
% Summarize the RGC layer properties in a string
%
% Syntax:
%   str = describe(obj)
%   str = @rgcLayer.describe;
%
% Description:
%    Provide a summary string of the RGC layer properties.
%
% Inputs:
%    obj - Object. A rgc layer object.
%
% Outputs:
%    str - String. The summary string.
%
% Optional key/value pairs:
%    None.
%

% History:
%    XX/XX/17  BW   ISETBIO Team, 2017
%    06/21/19  JNM  Documentation pass

%% Layer values
str = [];
txt = sprintf('Patch size:   \t %.1f %.1f um\n', obj.input.size * 1e6);
str = addText(str, txt);
txt = sprintf('Center:       \t   %.1f, %.1f mm\n', ...
    [obj.input.center] * 1e3);
str = addText(str, txt);
txt = sprintf('Time step:    \t %.1f ms\n', ...
    obj.input.input.integrationTime * 1e3);
str = addText(str, txt);
txt = sprintf('Duration:     \t %.1f ms\n', ...
    1e3 * obj.timeStep * (size(obj.mosaic{1}.responseLinear, 3)));
str = addText(str, txt);
txt = sprintf('N Trials:     \t %.0f \n', obj.nTrials);
str = addText(str, txt);
txt = sprintf('Eye side:     \t %s \n', obj.input.input.whichEye);
str = addText(str, txt);
txt = sprintf('Species:      \t %s \n', obj.input.input.species);
str = addText(str, txt);
str = addText(str, sprintf('---\n'));

%% Selected mosaic values
gdata = guidata(obj.fig);
nMosaic = get(gdata.listMosaics, 'Value');
[r, c, ~] = size(obj.mosaic{nMosaic}.cellLocation);
txt = sprintf('Row, Col:     \t %d, %d\n', r, c);
str = addText(str, txt);

end