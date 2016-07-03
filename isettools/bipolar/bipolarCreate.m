function obj = bipolarCreate(os, varargin)
% Create a bipolar object
%
%   bp = bipolarCreate(outerSegment);
%
% See @bipolar for more details on the bipolar class.
% 
% Inputs:  outer segment object
%
% Outputs:    the bipolar object
%
% Examples:
%   bp = bipolarCreate();
%   bp = bipolarCreate(os);
%
% JRG ISETBIO Team, Copyright, 5/2016

p = inputParser; 

addRequired(p, 'os');
addParameter(p, 'type', 'subunit')
% Parse and put results into structure p.
p.parse(os, varargin{:}); params = p.Results;
os = params.os;

%% Create the proper object
switch ieParamFormat(params.type)
    case {'linear'}
        % obj = bipolarLinear();
        warning('Not yet implemented');
    case {'subunit'}
        obj = bipolar(os);
    otherwise
        obj = bipolar(os);
end

end

