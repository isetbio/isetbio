function str = describe(obj)
% DESCRIBE - Summarize the bipolar layer properties in a string
%
%   @bipolarLayer.describe;
%   
% BW ISETBIO Team, 2017

%%
gdata = guidata(obj.figureHandle);
nMosaic = get(gdata.listMosaics,'Value');

%% General properties
str = [];
txt = sprintf('N Mosaics:    \t%d \n',length(obj.mosaic));
str = addText(str,txt);
txt = sprintf('Patch size:   \t%.1f %.1f um\n',obj.input.size*1e6);
str = addText(str,txt);
txt = sprintf('Center:       \t(%.1f,%.1f) mm\n',[obj.input.center]*1e3);
str = addText(str,txt);
txt = sprintf('Time step:    \t%.1f ms\n',obj.input.integrationTime*1e3);
str = addText(str,txt);
txt = sprintf('N Trials:     \t%.0f \n',obj.numberTrials);
str = addText(str,txt);
txt = sprintf('Eye side:     \t%s \n',obj.input.whichEye);
str = addText(str,txt);
txt = sprintf('Species:      \t%s \n',obj.input.species);
str = addText(str,txt);

% Specific to this mosaic
txt = sprintf('---\n');
str = addText(str,txt);
[r,c,~] = size(obj.mosaic{nMosaic}.responseCenter);
txt = sprintf('Samples:   \t%d %d \n',r,c);
str = addText(str,txt);

end