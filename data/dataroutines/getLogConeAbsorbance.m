function [log10absorbance,wave,params,comment] = coneAbsorbanceReadData(varargin)
%%coneAbsorbanceReadData  Return the log10 of the cone photopigment absorbance (aka optical density)
%
% Syntax:
%    log10absorbance = coneAbsorbanceReadData;
%
% Description:
%    Return the log10 of the cone photopigment absorbance (aka optical density).
%
%    Do not change the defaults for this routine without being aware that they are
%    used by other routines (e.g. photoPigment by default gets its absorbance via
%    this routine.)
%
%    It's not entirely clear that its best to have this returned as the log10 of
%    the quantity we ultimately use, but that is the format that Stockman stored
%    it on on his web page and there are advantages to staying with that.
%
% Input:
%    None.
%
% Output:
%     log10absorbance            Log10 of the cone absorbance, in the columns of a matrix.
%
%     wave                       Column vector of sample wavelengths, in nm.
%
%     params                     Structure of key/value pairs used to generate data.
%
%     comment                    A short comment describing the data, returned as a string.
%
% Optional key/value pairs:
%     'species'                  String specifying species
%                                  'human' (default) - Human L, M and S cones.
%                                  'monkey' - Monkey L, M and S cones.  Currently the same as human.
%
%     'coneAbsorbanceSource'     Source of the data
%                                  'StockmanSharpe' (default).
%                                   Values are StockmanSharpe estimates.  Valid when 'species' is 'human' or
%                                   'monkey'.  Data taken from Psychtoolbox mat file T_log10coneabsorbance_ss.
%                                   These in turn came from CVRL (http://cvrl.org).
%
%                                  The value for 'coneAbsorbanceSource' may be passed as a function handle, in
%                                  which case the passed function is called direclty with the key/value pairs passed to this
%                                  routine. 
%
%     'wave'                     Column vector of evenly spaced sample wavelengths in nm (default, (390:830)').
%
% See also: getRawData

% 08/10/17  dhb  Drafted.

%% Parse inputs
p = inputParser;
p.KeepUnmatched = true;
p.addParameter('species','human', @ischar);
p.addParameter('coneAbsorbanceSource','StockmanSharpe', @ischar);
p.addParameter('wave',(390:830)', @isnumeric);
p.parse(varargin{:});
params = p.Results;

%% Take care of case where a function handle is specified as source
%
% This allows for custom data to be defined by a user, via a function that
% could live outside of ISETBio.
if (isa(params.coneAbsorbanceSource,'function_handle'))
    [log10absorbance,wave,params,comment] = params.coneAbsorbanceSource(varargin{:});
    return;
end

%% Handle choices
switch (params.species)
    case {'human', 'monkey'}
        switch (params.coneAbsorbanceSource)
            case 'StockmanSharpe'
                % Load the absorbance from the PTB style T_ file, spline to desired wavelength
                % sampling, transpose to match ISETBio convention, and return.
                theData = getRawData('T_log10coneabsorbance_ss','datatype','ptbmatfileonpath');
                wavein = SToWls(theData.S_log10coneabsorbance_ss);
                wave = params.wave;
                log10absorbance = SplineCmf(wavein,theData.T_log10coneabsorbance_ss,wave)';
                comment = 'Log10 cone absorbance (aka optical density) estimated by Stockman/Sharpe, from CVRL via PTB';
            otherwise
                error('Unsupprted source specified');
        end
        
    otherwise
        error('Unsupported species specified');
end
        