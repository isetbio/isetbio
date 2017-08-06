function [uData,h] = plot(os, pType, varargin)
% Plot for the outersegment base class
%
% Some special cases of plot are managed in osLinear, osBioPhys.  General
% cases, however, are sent to this routine as
% plot@outersegment(obj,varargin).
%
% Inputs: 
%   os - outersegment object (any time)
% 
% pType Options:
%   isomerizations
%   current
%   all
%
% Outputs: 
%   uData - userdata that was plotted (uData = get(gca,'userdata'))
%   h is a handle to the plot window
%
% Properties that can be plotted:
%   None yet.
%
% Examples:
%
% JRG/BW (c) isetbio

%% Check for the number of arguments and create parser object.

p = inputParser;
p.addRequired('os');
p.addRequired('pType', @ischar);

p.parse(os, pType, varargin{:});
% No parameter values set yet

uData = []; h = [];

%% Choosing the plot type
switch ieParamFormat(pType)
           
    otherwise       
        error('Unknown outersegment base class plot type %s\n',pType);
end

end



