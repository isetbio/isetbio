function [ois, scene] = ...
    oisCreate(oisType, composition, modulation, varargin)
% Convenient script for common oi sequences
%
% Syntax
%   [ois, scene] = oisCreate(oisType, composition, modulation, [varargin])
%
% Description:
%    An oiSequence specifies certain simple retinal images that vary over
%    time, such as stimuli used in typical psychophysical experiments. This
%    function produces several common oiSequences.
%
%    The sequence is a mixture of a fixed OI and a modulated OI. The
%    mixture is determined by a time series of weights. The weights are
%    used for a mixture that is either an addition
%
%        oiFixed + w * oiModulated,
%
%    or a blend
%
%        w * oiFixed + (1 - w) * oiModulated
%
%    There are examples contained in the code. To access, type 'edit
%    oisCreate.m' into the Command Window.
%
% Inputs:
%    oisType         - String. One of 'vernier', 'harmonic', 'impulse'.
%    composition     - String. Either 'add' or 'blend'.
%    modulation      - Vector. A series of weights describing the add or
%                      blend composition levels.
%
% Outputs:
%    ois             - Struct. The created oi sequence structure.
%    scene           - Struct. A scene structure.
%
% Optional key/value pairs:
%    sampleTimes     - Vector. A vector of sample times. Default [].
%    testParameters  - Struct. A parameter structure for test targets.
%                      Default [].
%    sceneParameters - Struct. A scene parameter structure
%                      (e.g., fov, luminance). Default [].
%    oi              - Struct. An optical image structure. Default [].
%
% References:
%   ISETBIO wiki: https://github.com/isetbio/isetbio/wiki/OI-Sequences
%
% See Also:
%   t_oisCreate, oiSequence, sceneCreate, sceneHarmonic,
%   humanConeContrast, humanConeIsolating
%

% History:
%    xx/xx/16  BW   ISETBIO Team, 2016
%    03/19/18  jnm  Formatting
%    06/28/19  JNM  Formatting update

% Examples:
%{
    % Harmonic
    clear hparams
    hparams(2) = harmonicP;
    hparams(2).freq = 6;
    hparams(2).GaborFlag = .2;
    hparams(1) = hparams(2);
    hparams(1).contrast = 0;
    sparams.fov = 1;
    stimWeights = ieScale(fspecial('gaussian', [1, 50], 15), 0, 1);
    ois = oisCreate('harmonic', 'blend', stimWeights, ...
        'testParameters', hparams, 'sceneParameters', sparams);
    ois.visualize('movie illuminance');
%}
%{
    % Vernier
    clear vparams;
    vparams(2) = vernierP;
    vparams(2).name = 'offset';
    vparams(2).bgColor = 0;
    vparams(1) = vparams(2);
    vparams(1).barWidth = 0;
    vparams(1).bgColor = 0.5;
    vparams(1).name = 'uniform';
    sparams.fov = 1;
    stimWeights = ieScale(fspecial('gaussian', [1, 50], 15), 0, 1);
    [vernier, scenes] = oisCreate('vernier', 'add', stimWeights, ...
       'testParameters', vparams, 'sceneParameters', sparams);
    vernier.visualize('movie illuminance');

    ieAddObject(scenes{1});
    ieAddObject(scenes{2});
    sceneWindow;
%}
%{
    % Impulse (temporal)
    clear iparams
    sparams.fov = 1;
    sparams.meanluminance = 100;
    stimWeights = zeros(1, 50);
    stimWeights(2:4) = 1;
    impulse = oisCreate('impulse', 'add', stimWeights, ...
        'sceneParameters', sparams);
    impulse.visualize('movie illuminance');
%}

%% Inputs
p = inputParser;
p.KeepUnmatched = true;

% Required
validTypes = {'harmonic', 'vernier', 'impulse'};
validComp = {'add', 'blend'};
p.addRequired('oisType', @(x)(ismember(x, validTypes)));
p.addRequired('composition', @(x)(ismember(x, validComp)));
p.addRequired('modulation');

% Parameters
p.addParameter('sampleTimes', [], @isvector);
p.addParameter('testParameters', [], @isstruct);
p.addParameter('sceneParameters', [], @isstruct);
p.addParameter('oi', [], @isstruct);

p.parse(oisType, composition, modulation, varargin{:});

oi = p.Results.oi;
if isempty(oi)
    oi = oiCreate('wvf human');
    oi = oiSet(oi,'compute method','opticsotf');
    disp('Creating for opticsOTF method.')
end

oisType = ieParamFormat(oisType);

sampleTimes = p.Results.sampleTimes;
if isempty(sampleTimes)
    % By default, sample times are spaced 1 ms apart over the modulation.
    sampleTimes = 0.001 * ((1:length(modulation)) - 1);
end

% Test stimulus parameters (two structs needed for harmonic/vernier)
tparams = p.Results.testParameters;
% Scene parameters applied to both scenes (fov, luminance, etc.)
sparams = p.Results.sceneParameters;

%%
switch oisType
    case 'harmonic'
        if length(tparams) ~= 2
            error('Specify two harmonic param sets.');
        end
        scene = cell(1, 2);
        OIs   = cell(1, 2);

        % Color case requires specification of wave for SPDs
        if isfield(tparams(1), 'wave'), wave = tparams(1).wave;
        else, wave = 400:10:700;
        end

        for ii = 1:2
            scene{ii} = sceneCreate('harmonic', tparams(ii), wave);
            sname = sprintf('F %d C %.2f', tparams(ii).freq, ...
                tparams(ii).contrast);
            scene{ii} = sceneSet(scene{ii}, 'name', sname);
        end
        scene = localApplySparams(scene, sparams);

        for ii = 1:2
            OIs{ii} = oiCompute(oi, scene{ii}, 'pad value', 'mean');
        end
        ois = oiSequence(OIs{1}, OIs{2}, sampleTimes, modulation, ...
            'composition', composition);

    case 'vernier'
        if length(tparams) ~= 2
            error('Specify two vernier param sets.');
        end
        scene = cell(1, 2);
        OIs   = cell(1, 2);

        for ii = 1:2
            scene{ii} = sceneVernier('vernier', 'display', tparams(ii));
            scene{ii} = sceneSet(scene{ii}, 'name', tparams(ii).name);
        end
        scene = localApplySparams(scene, sparams);

        for ii = 1:2
            OIs{ii} = oiCompute(oi, scene{ii}, 'pad value', 'mean');
        end
        ois = oiSequence(OIs{1}, OIs{2}, sampleTimes, modulation, ...
            'composition', composition);

    case 'impulse'
        scene = cell(1, 2);
        OIs   = cell(1, 2);

        for ii = 1:2, scene{ii} = sceneCreate('uniform ee'); end
        scene = localApplySparams(scene, sparams);

        for ii = 1:2
            OIs{ii} = oiCompute(oi, scene{ii}, 'pad value', 'mean');
        end
        ois = oiSequence(OIs{1}, OIs{2}, sampleTimes, modulation, ...
            'composition', composition);

    otherwise
        error('Unknown type %s\n', oisType);
end

end

% -------------------------------------------------------------------------

function scenes = localApplySparams(scenes, sparams)
% Apply fields of sparams struct to every scene in the cell array.
if isempty(sparams), return; end
fields = fieldnames(sparams);
for ii = 1:numel(fields)
    for jj = 1:numel(scenes)
        scenes{jj} = sceneSet(scenes{jj}, fields{ii}, sparams.(fields{ii}));
    end
end
end
