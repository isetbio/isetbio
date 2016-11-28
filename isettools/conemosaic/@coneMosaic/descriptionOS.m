function str = descriptionOS(obj,varargin)
% Summarize properties of the outer segment
%
%   obj:  Cone mosaic object
% 
% Create a string to display in the outer segment panel of the cone mosaic
%
% BW, ISETBIO Team, 2016

os = obj.os;

% Which type of outer segment
str = sprintf('OS:     %s\n',class(os));

% Noise flag status
if os.noiseFlag, noiseFlag = 'on';
else noiseFlag = 'off';
end
txt = sprintf('noise:  %s',noiseFlag);
str = addText(str,txt);

end
