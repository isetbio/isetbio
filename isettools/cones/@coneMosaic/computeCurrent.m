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
p.addParameter('bgR',[],@isnumeric);
p.KeepUnmatched = true;
p.parse(varargin{:});
bgR = p.Results.bgR;

<<<<<<< HEAD
% EK: Check that absorption time series has been computed
if (isempty(obj.absorptions)  || size(obj.absorptions,3) == 1) && (isempty(p.Results.absorptionsInXWFormat))
    error('You must compute isomerizations (absorptions) time series prior to the current.');
=======
% Check that absorption time series has been computed
if (isempty(obj.absorptions)  || size(obj.absorptions, 3) == 1) && ...
        (isempty(p.Results.absorptionsInXWFormat))
    error(['You must compute isomerizations (absorptions) time series ' ...
        'prior to the current.']);
>>>>>>> 24af784f75526c07d43761aa0613a2984fc579f7
end

% This is the background absorption rate. We pass it in to 'warm up' the
% biophysical model to reach steady state faster. It is also used by the
% linear os model to obtain the needed filters.
<<<<<<< HEAD
=======
bgR = coneMeanIsomerizations(obj, 'absorptionsInXWFormat', ...
    p.Results.absorptionsInXWFormat);
>>>>>>> 24af784f75526c07d43761aa0613a2984fc579f7

% EK: if background absorption rate is already computed and defined as input
% variable, use that bgR. If not, compute it on the spot..
if isempty(bgR)
    bgR = coneMeanIsomerizations(obj, 'absorptionsInXWFormat', p.Results.absorptionsInXWFormat);
end
%% Call the appropriate outer segment photocurrent computation
<<<<<<< HEAD
if isa(obj.os,'osLinear')
    [obj.current, interpFilters, meanCur] = obj.os.osCompute(obj,'bgR',mean(bgR),varargin{:});
elseif isa(obj.os,'osBioPhys')
    obj.current = obj.os.osCompute(obj,'bgR',bgR); % EK: previously mean(bgR), now changed to bgR, since that is already an average, to warm up the biophys model.
=======
if isa(obj.os, 'osLinear')
    [obj.current, interpFilters, meanCur] = ...
        obj.os.osCompute(obj, 'bgR', mean(bgR), varargin{:});
elseif isa(obj.os, 'osBioPhys')
    obj.current = obj.os.osCompute(obj, 'bgR', mean(bgR), varargin{:});
>>>>>>> 24af784f75526c07d43761aa0613a2984fc579f7
    interpFilters = [];
    meanCur = [];
else
    error('Attempting to computer current with unsupported os class');
end

end