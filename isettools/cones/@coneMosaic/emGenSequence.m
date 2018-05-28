function [emPath, fixEMobj] = emGenSequence(obj, nEyeMovements, varargin)
% Generate sequence of fixational eye movements for a rect cone mosaic
%
% Syntax:
%	[emPath, fixEMojb] = emGenSequence(obj, nEyeMovements, [varargin])
%
% Description:
%  The eye movement samples are created at the same temporal sample
%  rate as the cone integration time. We only update the position at
%  the beginning of each integration time.
%
% Inputs:
%     obj               - rect cone mosaic object
%     nEyeMovements     - number of eye movements to generate
%
% Ouputs:
%     emPath            - eye positions in a tensor of size 
%                            nTrials x nEyeMovements x 2 matrix.  
%                         Units are cone positions
%     fixEMobj          - fixational eye movement object containing all the
%                         parameters
%
% Optional key/value pairs:
%    'em'               - Fixational eye movement structure, see
%                         fixationalEM for details; if not passed in
%                         the default fixationalEM is used.
%    'microsaccadetype' - One of these types of microsaccdade models
%        'none'  -  (default)
%        'stats based'
%        'heatmap/fixation based'
%    'rSeed'            - Random seed to be used
%    'nTrials'          - Multiple trial case, default = 1
%
%
% See Also:
%     fixationalEM
%

% History:
%    xx/xx/16  HJ/BW    ISETBIO Team, 2016
%    11/06/17  ncp      Added line to make drift magnitude independent of
%                       sample time.
%    11/06/17  dhb/npc  Added comments on microsaccade algorithm.
%    11/07/17  dhb      More cleaning and robustness.
%    02/26/18  jnm      Formatting, fix example
%    05/13/18  baw      Re-wrote for fixationalEM class.

% Examples:
%{
 cm = coneMosaic;
 cm.emGenSequence('help');
%}
%{
 scene = sceneCreate('mackay'); scene = sceneSet(scene,'fov',1);
 oi = oiCreate; oi = oiCompute(oi,scene);
 cm = coneMosaic; 
 cm.emGenSequence(50, 'nTrials', 1);
 cm.compute(oi);  cm.window;
%}
%{
 scene = sceneCreate('mackay'); scene = sceneSet(scene,'fov',1);
 oi = oiCreate; oi = oiCompute(oi,scene);
 cm = coneMosaic; 
 cm.emGenSequence(50, ...
        'nTrials', 1, ...
        'microsaccade type',...
        'heatmap/fixation based');
 cm.compute(oi);  cm.window;
%}

%% Help
if strcmp(nEyeMovements,'help')
 doc('coneMosaic.emGenSequence');
 return;
end

%% parse inputs
p = inputParser;
varargin = ieParamFormat(varargin);

p.addRequired('obj', @(x)(isa(x,'coneMosaic')));
p.addRequired('nEyeMovements', @isscalar);

% Either create the default, or the user creates it with special parameters
% and provides it
p.addParameter('em', fixationalEM, @(x)(isa(x,'fixationalEM')));

% User can set the microSaccadeType, but nothing else, about the
% fixationalEM.
validTypes = {'none','stats based','heatmap/fixation based'};
p.addParameter('microsaccadetype', 'none', @(x)(ismember(x,validTypes)));
p.addParameter('ntrials', 1, @isscalar);
p.addParameter('rseed', 1, @isscalar);
p.addParameter('computevelocity', [], @islogical);
p.addParameter('useparfor', false, @islogical);

% set parameters
p.parse(obj, nEyeMovements, varargin{:});

fixEMobj   = p.Results.em;
nTrials    = p.Results.ntrials;
randomSeed = p.Results.rseed;
if ~isempty(p.Results.rseed), rng(p.Results.rseed); end

microSaccadeType = p.Results.microsaccadetype;

%% Start the calculation

fixEMobj.microSaccadeType = microSaccadeType;
fixEMobj.randomSeed = randomSeed;
fixEMobj.computeForConeMosaic(obj,nEyeMovements, ...
    'nTrials',nTrials, ...
    'rSeed', randomSeed);

%% Set the cone eye movement positions variable

obj.emPositions = fixEMobj.emPos;
if nargout > 0, emPath = fixEMobj.emPos; end

end

