function [l,m,s] = coneTypeLocations(cmosaic,varargin)
% CONTYPELOCATIONS - Return row/col or indices of the three types of cones
%
%   [l,m,s] = coneTypeLocations(cmosaic, 'val', {'index' or 'rowcol'})
%
% Inputs
%   cmosaic - Cone Mosaic object
%  
% Parameter-Value inputs
%   'format' - {rowcol, index}  Return indices (default) or rowcol values
% 
% The l,m,s outputs are either 
%
%    an index to each cone type  (default), or
%    an N x 2 matrix of the row,col locations of the l, m, and s cone types
%    
% Example:
%    cm = coneMosaic;
%    [l,m,s] = coneTypeLocations(cm);
%    [l,m,s] = coneTypeLocations(cm,'format','rowcol');
%
% JRG/BW ISETBIO Team, 2016

%%
p = inputParser;
p.addRequired('cmosaic',@(x)(isa(cmosaic,'coneMosaic')));
p.addParameter('format','index',@ischar);
p.parse(cmosaic,varargin{:});

format = p.Results.format;

%% Compute

% Always compute indices
l  = find(cmosaic.pattern == 2);
m  = find(cmosaic.pattern == 3);
s  = find(cmosaic.pattern == 4);

switch format
    case 'index'
        % Done
    case 'rowcol'       
        % Get the (row,col) locations of each of the cone type
        [r,c] = ind2sub(size(cmosaic.pattern),l); 
        l = [r,c]; 
        [r,c] = ind2sub(size(cmosaic.pattern),m); 
        m = [r,c];
        [r,c] = ind2sub(size(cmosaic.pattern),s); 
        s = [r,c];
    otherwise
        error('Unknown return value type %s\n',format);
end

end
