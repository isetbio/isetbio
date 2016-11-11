function computeCurrent(obj, varargin)
% Convert absorptions to photocurrent with the os model.
% 
% obj:  A coneMosaic
% 
% HJ ISETBIO Team 2016

% parse inputs
p = inputParser;
p.KeepUnmatched = true;
% p.addParameter('bgR', 0, @isnumeric);

p.parse(varargin{:});
% bgR = p.Results.bgR;

% We should check that absorptions is not empty
if isempty(obj.absorptions)
    disp('Compute absorptions first.  No current comnputed');
end

% Pass in the cone mosaic object and any parsed arguments
obj.os.osCompute(obj, p.Results);

end