function dCell = cellDelete(c, dList)
% Delete some cells from a list of cells 
%
% Syntax:
%   dCell = cellDelete(c, dList);
%
% Description:
%    Delete a the cell entries in the integer list, dList. Return the
%    (shortened) cell array. 
%
% Inputs:
%    c     - The cell array
%    dList - The cell entries to remove
%
% Outputs:
%    dCell - The modified cell array
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    12/12/17  jnm  Formatting
%    01/19/18  jnm  Formatting update to match Wiki.

% Example:
%{
	a = {'a', 'b', 'c', 'd'};
    b = cellDelete(a, [1, 3]);
%}

if notDefined('c'), error('Cell required.'); end
if notDefined('dList'), error('Delete list required.'); end

if max(dList) > length(c) || min(dList < 1), error('Bad dList'); end

% there is actually a simple and much faster way to do that:
dCell = c;
dCell(dList)=[];
    
end