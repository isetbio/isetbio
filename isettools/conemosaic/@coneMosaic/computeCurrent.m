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

% We check that absorptions is not empty
if isempty(cMosaic.absorptions)
    disp('You must compute absorptions first.  Although we could do it here.');
end

% OLD:
%  obj.os.osCompute(obj.absorptions/obj.integrationTime, obj.pattern, 'append', false);

%% Make the call to the outer segment photocurrent computation

% This is the background rate we expect.  By passing it in, we warm up the
% biophysical model to reach steady state faster.  
bgR = coneMeanIsomerizations(cMosaic);

% Onward.
cMosaic.os.osCompute(cMosaic,'bgR',mean(bgR),varargin{:});

end