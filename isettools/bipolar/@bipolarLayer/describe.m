function str = describe(obj)
% DESCRIBE - Summarize the bipolar layer properties in a string
%
%   @bipolarLayer.describe;
%   
% BW ISETBIO Team, 2017

str = [];
txt = sprintf('N Mosaics:    \t%d \n',length(obj.mosaic));
str = addText(str,txt);
txt = sprintf('Patch size:   \t%.1f um\n',obj.size*1e6);
str = addText(str,txt);
txt = sprintf('Center:       \t(%.1f,%.1f) mm\n',[obj.center]*1e3);
str = addText(str,txt);
txt = sprintf('Time step:    \t%.1f ms\n',obj.timeStep*1e3);
str = addText(str,txt);
txt = sprintf('N Trials:     \t%.0f \n',obj.numberTrials);
str = addText(str,txt);
txt = sprintf('Eye side:     \t%s \n',obj.eyeSide);
str = addText(str,txt);
txt = sprintf('Species:      \t%s \n',obj.species);
str = addText(str,txt);
            
end