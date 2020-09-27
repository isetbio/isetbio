function [uData, h] = plot(os, pType, varargin)
% Plot for the outersegment base class
%
% Syntax:
%   [uData, h] = plot(os, pType, [varargin])
%
% Description:
%    Some special cases of plot are managed in osLinear, osBioPhys. General
%    cases, however, are sent to this routine as:
%
%         plot@outersegment(obj, varargin).
%
%    Properties that can be plotted:
%         None yet.
% Inputs: 
%	 os       - outersegment object (any time)
%    pType    - Options include: isomerizations, current, all.
%    varargin - (Optional) Additional arguments. **More information
%               required here**
%
% Outputs: 
%   uData     - Userdata that was plotted (uData = get(gca, 'userdata'))
%   h         - A handle to the plot window
%
% Optional key/value pairs:
%    None(?).
%

% History:
%    xx/xx/xx  JRG/BW  (c) isetbio
%    02/12/18  jnm     Formatting

%% Check for the number of arguments and create parser object.

p = inputParser;
p.addRequired('os');
p.addRequired('pType', @ischar);

p.parse(os, pType, varargin{:});
% No parameter values set yet

uData = [];
h = [];

%% Choosing the plot type
switch ieParamFormat(pType) 
    otherwise       
        error('Unknown outersegment base class plot type %s\n', pType);
end

end
