function merged = cellMerge(varargin)
%Merge a set of cell arrays into one long cell array
%
%   merged = cellMerge(varargin);
%
%Purpose:
%    Merge a list of cell arrays into one long cell array.
%
% Example:
%  a = {'a','b'}; b = {'c'}; cellMerge(a,b);
%  a = {'a','b'}; b = {'c'}; c = []; cellMerge(a,b,c)
%
% Copyright ImagEval Consultants, LLC, 2003.

warning('Deprecated. Use [cell1 cell2] instead.');
if isempty(varargin), merged = []; return; end

n = length(varargin);
for ii=1:n
    if ~isempty(varargin{ii}) && ~iscell(varargin{ii}), 
        error('All arguments must be cell arrays'); 
    end
    eachLength(ii) = length(varargin{ii});
end

totalLength = sum(eachLength);

this = 1;
for ii=1:n
    % Empty cell arrays are permitted
    if eachLength(ii) > 0
        for jj=1:eachLength(ii)
            merged{this} = varargin{ii}{jj};
            this = this+1;
        end
    end
end

end