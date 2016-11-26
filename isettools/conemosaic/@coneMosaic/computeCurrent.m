function interpFilters = computeCurrent(cMosaic, varargin)
% Convert absorptions to photocurrent using the os model.
% 
% Input:
%   cMosaic:  A coneMosaic
%
% Return
%   interpFilters - The linear filters are usually computed to a 1 ms time
%       base from the biophys model.  But they are interpolated to the time
%       samples of the cone mosaic.  The interpolated filters are provided
%       here.  To get the 1 ms values, use osLinear.lineFilters;
%
% HJ ISETBIO Team 2016

%% parse inputs
p = inputParser;
p.KeepUnmatched = true;
p.parse(varargin{:});

% Check that the absorptions have been computed
if isempty(cMosaic.absorptions)
    disp('You must compute isomerizations (absorptions) prior to the current.');
end

interpFilters = [];

%% Call the relevant outer segment photocurrent computation

% This is the background rate we expect.  By passing it in, we warm up the
% biophysical model to reach steady state faster.  
bgR = coneMeanIsomerizations(cMosaic);

% If osLinear class
if isa(cMosaic.os,'osLinear')
    [cMosaic.current, interpFilters] = cMosaic.os.osCompute(cMosaic,'bgR',mean(bgR),varargin{:});
else
    cMosaic.current = cMosaic.os.osCompute(cMosaic,'bgR',mean(bgR),varargin{:});
end

end