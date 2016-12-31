function [interpFilters, meanCur] = computeCurrent(cMosaic, varargin)
% Convert absorptions to photocurrent using the os model.
% 
% Input:
%   cMosaic:  A coneMosaic
%
% Return
%   interpFilters - The linear filters are computed to a 1 ms or less time
%       base from the biophys model.  They are interpolated to the time
%       samples of the cone mosaic.  The interpolated filters are provided
%       here.  To get the 1 ms values, use osLinear.lineFilters;
%   meanCur  - To do the full computation we need the mean current, too.
%              That is 
%
% HJ ISETBIO Team 2016

%% parse inputs
p = inputParser;
p.KeepUnmatched = true;
p.parse(varargin{:});

% Optional returns for osLinear case
interpFilters = [];
meanCur       = [];

% Check that absorption time series has been computed
if isempty(cMosaic.absorptions)  || size(cMosaic.absorptions,3) == 1
    disp('You must compute isomerizations (absorptions) time series prior to the current.');
    return;
end

%% Call the appropriate outer segment photocurrent computation

% This is the background absorption rate.  We pass it in to 'warm up' the
% biophysical model to reach steady state faster.
bgR = coneMeanIsomerizations(cMosaic);

if isa(cMosaic.os,'osLinear')
    [cMosaic.current, interpFilters, meanCur] = cMosaic.os.osCompute(cMosaic,'bgR',mean(bgR),varargin{:});
    % ieMovie(cMosaic.current);
else
    % Biophysical class
    cMosaic.current = cMosaic.os.osCompute(cMosaic,'bgR',mean(bgR),varargin{:});
end

end