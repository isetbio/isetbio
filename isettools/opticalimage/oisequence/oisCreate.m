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
    sparams.luminance = 100;
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
if isempty(oi), oi = oiCreate('wvf human'); end

oisType = ieParamFormat(oisType);

sampleTimes = p.Results.sampleTimes;  % User defined sample times
if isempty(sampleTimes)
    % By default, the sample times are spaced 1ms apart and last for
    % the whole temporal modulation.
    sampleTimes = 0.001 * ((1:length(modulation)) - 1);
end

% Many of the types have a params structure that we will pass along
% Test stimulus parameters (two structs needed)
tparams = p.Results.testParameters;
% Scene parameters, only one set.
sparams = p.Results.sceneParameters;

%%
switch oisType
    case 'harmonic'
        % oisCreate('harmonic', ...);
        if length(tparams) ~= 2
            error('Specify two harmonic param sets.');
        end
        scene = cell(1, 2);
        OIs = cell(1, 2);

        % Color case requires specification of wave for SPDs
        if isfield(tparams(1), 'wave'), wave = tparams(1).wave;
        else, wave = 400:10:700;
        end

        % Create basic harmonics
        for ii = 1:2
            scene{ii} = sceneCreate('harmonic', tparams(ii), wave);
            sname = sprintf('F %d C %.2f', tparams(ii).freq, ...
                tparams(ii).contrast);
            scene{ii} = sceneSet(scene{ii}, 'name', sname);
        end

        % Adjust both scenes based on sparams.
        fields = fieldnames(sparams);
        if ~isempty(fields)
            for ii = 1:length(fields)
                for jj = 1:2
                    val = sparams.(fields{ii});
                    scene{jj} = sceneSet(scene{jj}, fields{ii}, val);
                end
            end
        end
        % ieAddObject(scene{1});
        % ieAddObject(scene{2});
        % sceneWindow;

        % Compute optical images from the scene
        for ii = 1:2, OIs{ii} = oiCompute(oi, scene{ii}); end
        % ieAddObject(OIs{1});
        % ieAddObject(OIs{2});
        % oiWindow;

        %% Build the oiSequence
        % The weights define some amount of the constant background and
        % some amount of the line on the same constant background
        ois = oiSequence(OIs{1}, OIs{2}, sampleTimes, modulation, ...
            'composition', composition);
        % ois.visualize('movie illuminance');

    case 'vernier'
        % oisCreate('vernier', ...);
        if length(tparams) ~= 2
            error('Specify two vernier param sets.');
        end
        scene = cell(1, 2);
        OIs = cell(1, 2);

        % Create vernier stimulus and background
        for ii = 1:2
            scene{ii} = sceneCreate('vernier', 'display', tparams(ii));
            scene{ii} = sceneSet(scene{ii}, 'name', tparams(ii).name);
        end

        % Adjust both scenes based on sparams.
        fields = fieldnames(sparams);
        if ~isempty(fields)
            for ii = 1:length(fields)
                for jj = 1:2
                    val = sparams.(fields{ii});
                    scene{jj} = sceneSet(scene{jj}, fields{ii}, val);
                end
            end
        end
        % ieAddObject(scene{1});
        % ieAddObject(scene{2});
        % sceneWindow;

        % Compute optical images from the scene
        for ii = 1:2, OIs{ii} = oiCompute(oi, scene{ii}); end
        % ieAddObject(OIs{1});
        % ieAddObject(OIs{2});
        % oiWindow;

        ois = oiSequence(OIs{1}, OIs{2}, sampleTimes, modulation, ...
            'composition', composition);
        % ois.visualize('movie illuminance');

    case 'impulse'
        % oisCreate('impulse', 'add', weights, 'sparams', sparams);
        scene = cell(1, 2);
        OIs = cell(1, 2);

        for ii = 1:2, scene{ii} = sceneCreate('uniform ee'); end

        % Adjust both scenes based on sparams.
        if ~isempty(sparams)
            fields = fieldnames(sparams);
            if ~isempty(fields)
                for ii = 1:length(fields)
                    for jj = 1:2
                        val = sparams.(fields{ii});
                        scene{jj} = sceneSet(scene{jj}, fields{ii}, val);
                    end
                end
            end
        end
        % ieAddObject(scene{1});
        % ieAddObject(scene{2});
        % sceneWindow;

        % Compute optical images from the scene
        for ii = 1:2, OIs{ii} = oiCompute(oi, scene{ii}); end
        % ieAddObject(OIs{1});
        % ieAddObject(OIs{2});
        % oiWindow;

        ois = oiSequence(OIs{1}, OIs{2}, sampleTimes, modulation, ...
            'composition', composition);
        % ois.visualize('movie illuminance');

    otherwise
        error('Unknown type %s\n', oisType);
end
