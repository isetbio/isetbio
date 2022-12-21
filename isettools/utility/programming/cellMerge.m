function merged = cellMerge(varargin)
% Merge a set of cell arrays into one long cell array
%
% Synopsis
%   merged = cellMerge(varargin);
%
%Purpose:
%    Merge a list of cell arrays into one long cell array.
%
% Example:
%  a = {'a','b'}; b = {'c'}; cellMerge(a,b);
%  a = {'a','b'}; b = {'c'}; c = []; cellMerge(a,b,c)
%
% See also
%   cellDelete
%

% Examples:
%{
a = {'a','b'};
b = {'c'};
cellMerge(a,b)
%}

if isempty(varargin), merged = []; return; end

n = length(varargin);
eachLength = zeros(1,n);
for ii=1:n
    if ~isempty(varargin{ii}) && ~iscell(varargin{ii})
        error('All arguments must be cell arrays');
    end
    eachLength(ii) = length(varargin{ii});
end

% totalLength = sum(eachLength);

this = 1;
for ii=1:n
    % Empty cell arrays are permitted
    if eachLength(ii) > 0
        for jj=1:eachLength(ii)
            merged{this} = varargin{ii}{jj}; %#ok<AGROW> 
            this = this+1;
        end
    end
end

end


