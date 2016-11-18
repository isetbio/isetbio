function params = harmonicP(varargin)
% Default harmonic params
%
%   p = harmonicP;
%
% To specify values for the params use param/val pairs
%
%
% BW, ISETBIO Team, 2016

%% Parse arguments

p = inputParser;

p.addParameter('ang',0,@isnumeric);
p.addParameter('contrast',1,@isnumeric);
p.addParameter('freq',1,@isnumeric);
p.addParameter('ph',pi/2,@isnumeric);
p.addParameter('row',64,@isscalar);
p.addParameter('col',64,@isscalar);
p.addParameter('GaborFlag',0,@isscalar);

p.parse(varargin{:});

%% Assign parameters

params.ang       = p.Results.ang;
params.contrast  = p.Results.contrast;
params.freq      = p.Results.freq; 
params.ph        = p.Results.ph;
params.row       = p.Results.row; 
params.col       = p.Results.col; 
params.GaborFlag = p.Results.GaborFlag;

end
