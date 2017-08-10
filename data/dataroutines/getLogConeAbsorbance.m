function [log10absorbance,wave,comment] = getLogConeAbsorbance(varargin)
%%getLogConeAbsorbance  Return the log10 of the cone photopigment absorbance (aka optical density).
%
% Syntax:
%    log10absorbance = getLogConeAbsorbance;
%
% Description:
%    Return the log10 of the cone photopigment absorbance (aka optical density).
%
% Input:
%    None.
%
% Output:
%     log10absorbance      Log10 of the cone absorbance, in the columns of a matrix.
%     wave                 Column vector of sample wavelengths, in nm.
%     comment              A short comment describing the data, returned as a string.
%
% Optional key/value pairs:
%     'species'            String specifying species
%                            'human' (default). Human L, M and S cones.
%                            'monkey' - Monkey L, M and S cones.  Currently the same as human.
%     'source'             Where the data are taken from
%                            'ptb' (default). Values from Psychtoolbox data (these in turn come from CVRL.)
%     'wave'               Column vector of evenly spaced sample wavelengths in nm (default, (390:830)').
%
% See also: getRawData

% 08/10/17  dhb  Drafted.

%% Parse inputs
p = inputParser;
p.addParameter('species','human', @ischar);
p.addParameter('source','ptb', @ischar);
p.addParameter('wave',(390:830)', @isnumeric);
p.parse(varargin{:});

%% Handle choices
switch (p.Results.species)
    case {'human', 'monkey'}
        switch (p.Results.source)
            case 'ptb'
                % Load the absorbance from the PTB style T_ file, spline to desired wavelength
                % sampling, transpose to match ISETBio convention, and return.
                theData = getRawData('T_log10coneabsorbance_ss','datatype','matfileonpath');
                wavein = SToWls(theData.S_log10coneabsorbance_ss);
                wave = p.Results.wave;
                log10absorbance = SplineCmf(wavein,theData.T_log10coneabsorbance_ss,wave)';
                comment = 'Log10 cone absorbance (aka optical density) estimated by Stockman/Sharpe, from CVRL via PTB';
            otherwise
                error('Unsupprted source specified');
        end
        
    otherwise
        error('Unsupported species specified');
end
        