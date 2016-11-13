function computeCurrent(cMosaic, varargin)
% Convert absorptions to photocurrent with the os model.
% 
% cMosaic:  A coneMosaic
% 
% HJ ISETBIO Team 2016

% parse inputs
p = inputParser;
p.KeepUnmatched = true;
p.parse(varargin{:});

% We should check that absorptions is not empty
if isempty(cMosaic.absorptions)
    disp('You must compute absorptions first.  Although we could do it here.');
end

% We should change this to pass in the cone mosaic object and any parsed
% arguments, rather than these parameters.
% obj.os.osCompute(obj.absorptions/obj.integrationTime, obj.pattern, 'append', false);
cMosaic.os.osCompute(cMosaic,varargin{:});

end