function parser = ieInputParser()
%% Get an inputParser configured with ISETBIO conventions
%
% parser = ieInputParser() returns an instance of the built-in inputParser
% class, with its properties set to ISETBIO defaults.
%
% See also inputParser
%
% parser = ieInputParser()
%
% Copyright (c) 2015 ISETBIO Team

parser = inputParser();

parser.CaseSensitive = false;
parser.KeepUnmatched = true;
parser.StructExpand = false;

% PartialMatching is not an option in Matlab < ~8.2/R2013b
%   That's OK.  We just want to turn it off when we have the choice.
if ~verLessThan('matlab', '8.2')
    parser.PartialMatching = false;
end

end
