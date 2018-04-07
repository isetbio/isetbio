function parser = ieInputParser()
% Get an inputParser configured with ISETBIO conventions
%
% Syntax:
%	parser = ieInputParser()
%
% Description:
%    returns an instance of the built-in inputParser class, with its
%    properties set to ISETBIO defaults.
%
% Inputs:
%    None.
%
% Outputs:
%    parser - The input parser object
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    inputParser
%

% History:
%    xx/xx/15       Copyright (c) 2015 ISETBIO Team
%    12/15/17  jnm  Formatting
%    01/19/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    parser = ieInputParser()
%}

parser = inputParser();

parser.CaseSensitive = false;
parser.KeepUnmatched = true;
parser.StructExpand = false;

% PartialMatching is not an option in Matlab < ~8.2/R2013b
% [Note: XXX - That's OK. We just want it off when we have the choice.]
if ~verLessThan('matlab', '8.2'), parser.PartialMatching = false; end

end
