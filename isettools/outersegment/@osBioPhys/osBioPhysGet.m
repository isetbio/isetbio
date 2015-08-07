function val = osBioPhysGet(obj, varargin)
% osBioPhysGet: a method of @osBioPhys that sets isetbio outersegment object 
% parameters using the input parser structure.
% 
% 
% 
% 
% 8/2015 JRG NC DHB

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Check for the number of arguments and create parser object.
%
% Here we make matching of key names case-insensitive, and cause
% any errors thrown in the parsing to report that they belong to this
% function, not to the parser object itself.
error(nargchk(0, Inf, nargin));
p = inputParser; p.CaseSensitive = false; p.FunctionName = mfilename;

% A set function always needs a string to say what is being set, and a
% value.  So we make those required arguments.
allowableFieldsToSet = {'noiseflag','conecurrentsignal'};
p.addRequired('what',@(x) any(validatestring(x,allowableFieldsToSet)));
% p.addRequired('value');

% Define what units are allowable
allowableUnitStrings = {'a', 'ma', 'ua', 'na', 'pa'}; % amps to picoamps

% Set up key value pairs
%
% Defaults are given here, as is a type check function
p.addParameter('units','pa',@(x) any(validatestring(x,allowableUnitStrings)));
% p.addParameter('sconefilter',0,@isnumeric);
% p.addParameter('mconefilter',0,@isnumeric);
% p.addParameter('lconefilter',0,@isnumeric);

% Do the parsing and put the results into an easy to access structure
p.parse(varargin{:}); params = p.Results;

% Dump out the parameters.  This wouldn't be here in real code, it's just
% so you can see what happens when you call this with various arguments.
% params

% And now there would be a nice big switch on params.what that would direct
% traffic to whatever needed doing.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

switch lower(params.what);  % Lower case and remove spaces

    case {'noiseflag'}
        
        val = obj.noiseFlag;
        
    case{'conecurrentsignal'}
        val = obj.ConeCurrentSignal;
        
end

