function [l,m,s] = coneTypeLocations(cmosaic,varargin)
% Return the row/col indices of the three types of cones
%
%   [l,m,s] = coneTypeLocations(cmosaic, 'val', {'index' or 'rowcol'})
%
% The l,m,s outputs are either 
%
%    an index to each cone type  (default), or
%    an N x 2 matrix of the row,col locations of the l, m, and s cone types
%    
% Example:
%    cm = coneMosaic;
%    [l,m,s] = coneTypeLocations(cm);
%
% JRG/BW ISETBIO Team, 2016

%%
p = inputParser;
p.addRequired('cmosaic',@(x)(isa(cmosaic,'coneMosaic')));
p.addParameter('val','index',@ischar);
p.parse(cmosaic,varargin{:});

val = p.Results.val;

%% Compute

% Always compute indices
l  = find(cmosaic.pattern == 2);
m  = find(cmosaic.pattern == 3);
s  = find(cmosaic.pattern == 4);

switch val
    case 'index'
        % Done
    case 'rowcol'       
        % Get the (row,col) locations of each of the cone type
        [r,c] = ind2sub(size(cmosaic.pattern),l); clear l m s; l = [r,c]; clear r c
        [r,c] = ind2sub(size(cmosaic.pattern),m); clear l m s; m = [r,c]; clear r c
        [r,c] = ind2sub(size(cmosaic.pattern),s); clear l m s; s = [r,c];
    otherwise
        error('Unknown return value type %s\n',val);
end

end
