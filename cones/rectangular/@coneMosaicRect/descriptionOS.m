function str = descriptionOS(obj)
% Summarize properties of the outer segment
%
% Syntax:
%   str = DESCRIPTIONSOS(obj)
%
% Description:
%    Create a string to display in the outer segment panel of the cone
%    mosaic, summarizing the properties.
%
%   Inputs:
%   obj - cone mosaic object
%
%   Outputs:
%   str - string describing some things about the mosaic's outer segment.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/16  BW   ISETBIO Team, 2016
%    02/19/18  jnm  Formatting

os = obj.os;

% Which type of outer segment
str = sprintf('OS:     %s\n', class(os));

% Noise flag status (random, frozen, none);
txt = sprintf('noise:  %s', os.noiseFlag);
str = addText(str, txt);

end
