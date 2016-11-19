function computeCurrent(cMosaic, varargin)
% Convert absorptions to photocurrent using the os model.
% 
%   cMosaic:  A coneMosaic
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

%% Call the relevant outer segment photocurrent computation

% This is the background rate we expect.  By passing it in, we warm up the
% biophysical model to reach steady state faster.  
bgR = coneMeanIsomerizations(cMosaic);

% Onward.
cMosaic.current = cMosaic.os.osCompute(cMosaic,'bgR',mean(bgR),varargin{:});

end