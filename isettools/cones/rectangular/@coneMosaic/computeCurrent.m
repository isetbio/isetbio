function [interpFilters, meanCur] = computeCurrent(obj, varargin)
% Convert absorptions to photocurrent using the os model.
%
% Syntax:
%   [interpFilters, meanCur] = computeCurrent(obj, [varargin]) 
%
% Input:
%    obj           - A coneMosaic object.
%
% Output:
%    interpFilters - The linear filters are computed from the biophys
%                    model. They are interpolated to the time samples of
%                    the cone mosaic. The interpolated filters are provided
%                    here. To get the 1 ms values, use
%                    osLinear.linearFilters.
%    meanCur       - Sometimes we need the mean current as well.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    CONEMOSAIC, COMPUTE, COMPUTEFOROISEQUENCE.
%

% History:
%    xx/xx/16  HJ   ISETBIO Team 2016
%    02/19/18  jnm  Formatting

%% parse inputs
p = inputParser;
p.addParameter('absorptionsInXWFormat', [], @isnumeric);
p.KeepUnmatched = true;
p.parse(varargin{:});

% Check that absorption time series has been computed
if (isempty(obj.absorptions)  || size(obj.absorptions, 3) == 1) && ...
        (isempty(p.Results.absorptionsInXWFormat))
    error(['You must compute isomerizations (absorptions) time series ' ...
        'prior to the current.']);
end

% This is the background absorption rate. We pass it in to 'warm up' the
% biophysical model to reach steady state faster. It is also used by the
% linear os model to obtain the needed filters.
bgR = coneMeanIsomerizations(obj, 'absorptionsInXWFormat', ...
    p.Results.absorptionsInXWFormat);

%% Call the appropriate outer segment photocurrent computation
if isa(obj.os, 'osLinear')
    [obj.current, interpFilters, meanCur] = ...
        obj.os.osCompute(obj, 'bgR', mean(bgR), varargin{:});
elseif isa(obj.os, 'osBioPhys')
    obj.current = obj.os.osCompute(obj, 'bgR', mean(bgR), varargin{:});
    interpFilters = [];
    meanCur = [];
else
    error('Attempting to computer current with unsupported os class');
end

end