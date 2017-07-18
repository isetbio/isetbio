function str = describe(obj)
% DESCRIBE - Summarize the RGC layer properties in a string
%
%   @rgcLayer.describe;
%   
% BW ISETBIO Team, 2017

str = [];
txt = sprintf('Patch size:   \t %.1f %.1f um\n',obj.input.size*1e6);
str = addText(str,txt);
txt = sprintf('Center:       \t   %.1f,%.1f mm\n',[obj.input.center]*1e3);
str = addText(str,txt);
txt = sprintf('Time step:    \t %.1f ms\n',obj.input.input.integrationTime*1e3);
str = addText(str,txt);
txt = sprintf('Duration:     \t %.1f ms\n',1e3*obj.timeStep*(size(obj.mosaic{1}.responseLinear,3)));
str = addText(str,txt);
txt = sprintf('N Trials:     \t %.0f \n',obj.nTrials);
str = addText(str,txt);
txt = sprintf('Eye side:     \t %s \n',obj.input.input.whichEye);
str = addText(str,txt);
txt = sprintf('Species:      \t %s \n',obj.input.input.species);
str = addText(str,txt);
            
end