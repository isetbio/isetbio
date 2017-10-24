function str = describe(obj)
% Syntax:
%
%   bipolar.describe;
%
% Dexcription:
%    Summarize the bipolar mosaic properties in a string
%
% Notes:
% * @JRG:  How will we handle multiple bipolar mosaics?
%

% History:
% BW ISETBIO Team, 2017
%
%    10/18/17  jnm  Comments & formatting

%% Summarize properties
str = sprintf('Cell type:     \t%s\n',obj.cellType);
txt = sprintf('Patch size:    \t%.1f um\n',obj.patchSize*1e6);
str = addText(str,txt);
txt = sprintf('Center size:     \t%.1f um\n',obj.sRFcenter(1));
str = addText(str,txt);
txt = sprintf('Surround size: \t%.1f um\n',obj.sRFsurround(1));
str = addText(str,txt);
txt = sprintf('Time step:       \t%.1f ms\n',obj.timeStep*1e3);
str = addText(str,txt);
txt = sprintf('Duration:        \t%.1f ms\n',obj.get('duration')*1e3);
str = addText(str,txt);
            
end
