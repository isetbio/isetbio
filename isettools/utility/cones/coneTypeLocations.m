function [l, m, s] = coneTypeLocations(cmosaic, varargin)
% Return row/col or indices of the three types of cones
%
% Syntax:
%   [l, m, s] = coneTypeLocations(cmosaic, 'val', [{'index' or 'rowcol'}])
%
% Description:
%    A function designed to return the row and column or indices of the
%    three types of cones (l, m, and s)
%
%    The l, m, s outputs are either 
%     - an index to each cone type  (default), or
%     - an Nx2 matrix of the row, col locations of the l, m, & s cone types
%
%    There are examples contained in the code. To access, type 'edit
%    coneTypeLocations.m' into the Command Window.
%
% Input:
%    cmosaic  - Cone Mosaic object
%
% Outputs:
%    l        - Pertinent information for L cones
%    m        - Pertinent information for M cones
%    s        - Pertinent information for S cones
%
% Optional key/value pairs:
%    'format' - {rowcol, index}  Return indices (default) or rowcol values
% 

% History:
%    xx/xx/16  JRG/BW  ISETBIO Team, 2016
%    02/09/18  jnm     Formatting

% Examples:
%{
   cm = coneMosaic;
   [l, m, s] = coneTypeLocations(cm);
   [l, m, s] = coneTypeLocations(cm, 'format', 'rowcol');
%}

%%
p = inputParser;
p.addRequired('cmosaic', @(x)(isa(cmosaic, 'coneMosaic')));
p.addParameter('format', 'index', @ischar);
p.parse(cmosaic, varargin{:});

format = p.Results.format;

%% Compute
% Always compute indices
l = find(cmosaic.pattern == 2);
m = find(cmosaic.pattern == 3);
s = find(cmosaic.pattern == 4);

switch format
    case 'index'
        % Done
    case 'rowcol'       
        % Get the (row, col) locations of each of the cone type
        [r, c] = ind2sub(size(cmosaic.pattern), l); 
        l = [r, c]; 
        [r, c] = ind2sub(size(cmosaic.pattern), m); 
        m = [r, c];
        [r, c] = ind2sub(size(cmosaic.pattern), s); 
        s = [r, c];
    otherwise
        error('Unknown return value type %s\n', format);
end

end
