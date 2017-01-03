function [interpFilters, meanCur] = computeCurrent(cMosaic, varargin)
%COMPUTECURRENT Convert absorptions to photocurrent using the os model.
%    [interpFilters, meanCur] = COMPUTECURRENT(cMosaic, varargin) 
%
% Input:
%   cMosaic:  A coneMosaic object.
%
% Output:
%   interpFilters - The linear filters are computed from the biophys model.
%       They are interpolated to the time samples of the cone mosaic.  The
%       interpolated filters are provided here.  To get the 1 ms values,
%       use osLinear.linearFilters.
%   meanCur - Sometimes we need the mean current as well.
%
% See also coneMosaic

% HJ ISETBIO Team 2016

%% parse inputs
p = inputParser;
p.KeepUnmatched = true;
p.parse(varargin{:});

% Check that absorption time series has been computed
if isempty(cMosaic.absorptions)  || size(cMosaic.absorptions,3) == 1
    error('You must compute isomerizations (absorptions) time series prior to the current.');
end

% This is the background absorption rate.  We pass it in to 'warm up' the
% biophysical model to reach steady state faster.  It is also used by the
% linear os model to obtain the needed filters.
bgR = coneMeanIsomerizations(cMosaic);

%% Call the appropriate outer segment photocurrent computation
if isa(cMosaic.os,'osLinear')
    [cMosaic.current, interpFilters, meanCur] = cMosaic.os.osCompute(cMosaic,'bgR',mean(bgR),varargin{:});
elseif isa(cMosaic.os,'osBioPhys')
    cMosaic.current = cMosaic.os.osCompute(cMosaic,'bgR',mean(bgR),varargin{:});
    interpFilters = [];
    meanCur       = [];
else
    error('Attempting to computer current with unsupported os class');
end

end