function str = descriptionOS(obj,varargin)
%DESCRIPTIONSOS  Summarize properties of the outer segment
%   str = DESCRIPTIONSOS(obj,varargin)
%
%   Create a string to display in the outer segment panel of the cone mosaic
%
%   Inputs:
%   obj     - cone mosaic object
%
%   Outputs:
%   str     - string describing some things about the mosaic.
%
%   Optional parameter name/value pairs chosen from the following:
%
%   None currently.

% BW, ISETBIO Team, 2016

os = obj.os;

% Which type of outer segment
str = sprintf('OS:     %s\n',class(os));

% Noise flag status (random, frozen, none);
txt = sprintf('noise:  %s',os.noiseFlag);
str = addText(str,txt);

end
