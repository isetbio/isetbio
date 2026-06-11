function exampleInputParser(varargin)
% Example parser using MATLAB to interpret key/value pairs.
%
% Syntax:
%   exampleInputParser(varargin)
%
% Description:
%    Example function that shows how to use Matlab's parser object to set
%    up a parameters structure using key/value pairs. This is set up to
%    show how we might do this for the set method for an object called
%    example.
%
%    The calling syntax is 
%        exampleInputParser('WhatToSet', valueToSet, ...)
%    where ... is set of key value pairs.
%
%    Supported key/value pairs
%       'units', 'unitstring'  - Set the units parameter to the string in
%                                unitstring (default, 'defaultunits').
%       'whatsi', whatsievalue - Set the whatsi parameter to the value in
%                                whatsivalue(default, 10)
%       'whosi', 'whosiname'   - Set the whosi parameter to the string in
%                                whosiname (default, 'joe')
%
%    The code illusrates how to check that the passed 'what' string
%    corresponds to an allowable value.
%
%    The code below also shows how to check whether passed 'unitstring' is
%    within a predefined set of allowable unit strings. You can get very
%    sophisticated about argument checking if you like. See the Matlab
%    documentation on inputParser, vaidatestring, and validateattributes.
%
%    With respect to specifying units, the parser by default returns the
%    string 'defaultunits'. The set of any particular field can then check
%    this and use the default units for that field if this is the string.
%    If some other units are explicitly specified, then the set code would
%    make sure that the units were allowable for its field and do the
%    conversion to the default units. Since a single set function may set
%    some parameters that are length and others that are time, etc. the
%    parser should probably do only a very broad check of allowable units,
%    or perhaps none at all.
%
%    In some early iset/isetbio functions, we allowed synonyms for some
%    keys, and also allowed keys to have whitespace. I am not sure we want
%    to continue with that, I prefer to just have one key for each
%    attribute, and I find it harder to read the high level code when I
%    don't know whether two different sets are really doing the same thing
%    or are doing different things.
% 
%    But if you wanted to do that, you could do so by pre-processing
%    varargin to map from the various synonyms onto the standard form
%    understood by the parser. And then the set code would just know about
%    the standard names of the resulting parameters fields.
%
%    It isn't clear how often we'll want key/value options other than units
%    for sets and gets, but in any case this routine shows how to do it.
% 
%    The input parsing for a get function would be similar, except that no
%    value would be passed.
%
% Inputs:
%    None required.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    Several key/value pair types are supported for this function. They
%    include the following:
%       'units', 'unitstring'  - Set the units parameter to the string in
%                                unitstring (default, 'defaultunits').
%       'whatsi', whatsievalue - Set the whatsi parameter to the value in
%                                whatsivalue(default, 10)
%       'whosi', 'whosiname'   - Set the whosi parameter to the string in
%                                whosiname (default, 'joe')
%

% History:
%    07/20/15  dhb  Wrote it. Sure beats working.
%    04/18/18  jnm  Formatting

% Examples:
%{
	exampleInputParser('oslength', 10);
	exampleInputParser('regenerationtime', 150, 'units', 'usec', ...
        'whatsi', 11, 'whosi', 'ellen');
%}

% Check for the number of arguments and create parser object.
%
% Here we make matching of key names case-insensitive, and cause
% any errors thrown in the parsing to report that they belong to this
% function, not to the parser object itself.
error(nargchk(0, Inf, nargin));
p = inputParser;
p.CaseSensitive = false;
p.FunctionName = mfilename;

% A set function always needs a string to say what is being set, and a
% value. So we make those required arguments.
allowableFieldsToSet = {'oslength' 'regenerationtime', 'somethinggain'};
p.addRequired('what', @(x) any(validatestring(x, allowableFieldsToSet)));
p.addRequired('value');

% Define what units are allowable
allowableUnitStrings = {'defaultunits', 'nm', 'um', 'mm', 'cm', 'm', ...
    'usec', 'msec', 'sec'};

% Set up key value pairs
%
% Defaults are given here, as is a type check function
p.addParameter('units', 'defaultunits', ...
    @(x) any(validatestring(x, allowableUnitStrings)));
p.addParameter('whatsi', 10, @isnumeric);
p.addParameter('whosi', 'joe', @ischar);

% Do the parsing and put the results into an easy to access structure
p.parse(varargin{:});
params = p.Results;

% Dump out the parameters. This wouldn't be here in real code, it's just
% so you can see what happens when you call this with various arguments.
params

% And now there would be a nice big switch on params.what that would direct
% traffic to whatever needed doing.

end
