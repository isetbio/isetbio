function str= describe(obj)
% Describe the RGC mosaic properties
%
% Prints the relevant text to a string, which is used in the display
% window. 
%
% BW, ISETBIO Team, 2017

parent = obj.Parent;  % Used for size and trials.  Needs help.

% Cell properties
str = sprintf('Cell type: %s\n',obj.cellType);
txt = sprintf('Model: %s\n',class(obj));
str = addText(str,txt);

% Mosaic properties
txt = sprintf('N Trials %d\n',parent.numberTrials);
str = addText(str,txt);
txt = sprintf('Patch size %d (um)\n',1e6*parent.size);
str = addText(str,txt);

% Spatial temporal properties
txt = sprintf('Row,Col,Time: %d, %d, %d\n',...
    size(obj.cellLocation),size(obj.responseLinear,3));
str = addText(str,txt);

txt = sprintf('Duration: %.1f ms\n',obj.dt*size(obj.responseLinear,3));
str = addText(str,txt);

end