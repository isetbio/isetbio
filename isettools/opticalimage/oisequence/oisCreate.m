function [ois, scene] = oisCreate(oisType,composition, modulation, varargin)
% OISCREATE - oi sequence creation
%    An oiSequence specifies certain simple retinal images that vary over
%    time, such as stimuli used in typical psychophysical experiments.
%
%    [ois, scenes] = OISCREATE(oisType,composition,modulation,'PARAM1',val ...)
%
%  Required parameters
%   'oisType'      - One of 'vernier','harmonic','impulse'.
%   'composition'  - 'add' or 'blend'
%   'modulation'   - Series of weights describing the add or blend
%
%  Optional parameter/val types chosen from the following 
%    'testParameters'   Parameters for the test targets 
%    'sceneParameters'  General scene parameters (e.g., fov, luminance)
%    
%    The sequence is a mixture of a fixed OI and a modulated OI. The
%    mixture is determined by a time series of weights.  The weights are
%    used for a mixture that is either an addition
%
%        oiFixed + w*oiModulated, 
%
%    or a blend
%
%        w*oiFixed + (1-w)*oiModulated 
%
%  Harmonics
%   clear hparams
%   hparams(2) = harmonicP; hparams(2).freq = 6; hparams(2).GaborFlag = .2; 
%   hparams(1) = hparams(2); hparams(1).contrast = 0; sparams.fov = 1; 
%   stimWeights = ieScale(fspecial('gaussian',[1,50],15),0,1);
%   ois = oisCreate('harmonic','blend',stimWeights, ...
%                   'testParameters',hparams, ...
%                   'ssceneParams',sparams);
%   ois.visualize;
%
%  Vernier
%   clear vparams; vparams(2) = vernierP; 
%   vparams(2).name = 'offset'; vparams(2).bgColor = 0; vparams(1) = vparams(2); 
%   vparams(1).barWidth = 0; vparams(1).bgColor = 0.5; vparams(1).name = 'uniform';
%   stimWeights = ieScale(fspecial('gaussian',[1,50],15),0,1);
%   [vernier, scenes] = oisCreate('vernier','add', stimWeights,...
%                                 'testParameters',vparams,...
%                                 'sceneParameters',sparams);
%   vernier.visualize;
%   ieAddObject(scenes{1}); ieAddObject(scenes{2}); sceneWindow;
%
% Impulse 
%   clear iparams
%   sparams.fov = 1; sparams.luminance = 100;
%   stimWeights = zeros(1,50); stimWeights(2:4) = 1;
%   impulse = oisCreate('impulse','add', stimWeights,...
%                       'sceneParameters',sparams);
%   impulse.visualize;
%
% See also SCENECREATE
%
% ISETBIO wiki: <a href="matlab:
%    web('https://github.com/isetbio/isetbio/wiki/OI-Sequences','-browser')">Optical image sequences</a>.
%
% BW ISETBIO Team, 2016

%% Inputs
p = inputParser;
p.KeepUnmatched = true;
p.addRequired('oisType',@ischar)
p.addRequired('composition',@ischar);
p.addRequired('modulation');

% Parameters that can be passed
p.addParameter('sampleTimes',[],@isvector);
p.addParameter('testParameters',[],@isstruct);
p.addParameter('sceneParameters',[],@isstruct);
p.addParameter('oi',[],@isstruct);

p.parse(oisType,composition,modulation,varargin{:});

oi = p.Results.oi;
if isempty(oi), oi = oiCreate('wvf human'); end

oisType = ieParamFormat(oisType);

sampleTimes = p.Results.sampleTimes;
if isempty(sampleTimes)
    sampleTimes = 0.001*((1:length(modulation))-1);
end

% Many of the types have a params structure that we will pass along
tparams = p.Results.testParameters;   % Test stimulus parameters (two structs needed)
sparams = p.Results.sceneParameters;   % Scene parameters, only one set.

%%
switch oisType
    case 'harmonic'
        % oisCreate('harmonic', ...); % See examples
        if length(tparams) ~= 2, error('Specify two harmonic param sets.'); end
        scene = cell(1,2);
        OIs = cell(1, 2);

        % Create basic harmonics
        for ii=1:2
            scene{ii} = sceneCreate('harmonic',tparams(ii));
            sname = sprintf('F %d C %.2f', tparams(ii).freq, tparams(ii).contrast);
            scene{ii} = sceneSet(scene{ii},'name',sname);
        end
        
        % Adjust both scenes based on sparams.
        fields = fieldnames(sparams);
        if ~isempty(fields)
            for ii=1:length(fields)
                for jj=1:2
                    val = sparams.(fields{ii});
                    scene{jj} = sceneSet(scene{jj}, fields{ii},val);
                end
            end
        end
        % ieAddObject(scene{1}); ieAddObject(scene{2}); sceneWindow;

        % Compute optical images from the scene
        for ii = 1:2
            OIs{ii} = oiCompute(oi,scene{ii});
        end
        % ieAddObject(OIs{1}); ieAddObject(OIs{2}); oiWindow;

        %% Build the oiSequence
        
        % The weights define some amount of the constant background and some amount
        % of the line on the same constant background
        ois = oiSequence(OIs{1}, OIs{2}, sampleTimes, modulation, ...
            'composition', composition);
        % ois.visualize;
    case 'vernier'
        % oisCreate('vernier', ...);   % See examples
        if length(tparams) ~= 2, error('Specify two vernier param sets.'); end
        scene = cell(1,2);
        OIs = cell(1, 2);
        
        % Create vernier stimulus and background
        for ii=1:2
            scene{ii} = sceneCreate('vernier', 'display', tparams(ii));
            scene{ii} = sceneSet(scene{ii},'name',tparams(ii).name);
        end
        
        % Adjust both scenes based on sparams.
        fields = fieldnames(sparams);
        if ~isempty(fields)
            for ii=1:length(fields)
                for jj=1:2
                    val = sparams.(fields{ii});
                    scene{jj} = sceneSet(scene{jj}, fields{ii},val);
                end
            end
        end
        %ieAddObject(scene{1}); ieAddObject(scene{2}); sceneWindow;

        % Compute optical images from the scene
        for ii = 1:2
            OIs{ii} = oiCompute(oi,scene{ii});
        end
        % ieAddObject(OIs{1}); ieAddObject(OIs{2}); oiWindow;

        ois = oiSequence(OIs{1}, OIs{2}, sampleTimes, modulation, ...
            'composition', composition);
        % ois.visualize;
        
    case 'impulse'
        % oisCreate('impulse', 'add', weights,'sparams',sparams); % See examples
        scene = cell(1,2);
        OIs = cell(1, 2);
        
        for ii=1:2
            scene{ii} = sceneCreate('uniform ee');
        end
        
        % Adjust both scenes based on sparams.
        if ~isempty(sparams)
            fields = fieldnames(sparams);
            if ~isempty(fields)
                for ii=1:length(fields)
                    for jj=1:2
                        val = sparams.(fields{ii});
                        scene{jj} = sceneSet(scene{jj}, fields{ii}, val);
                    end
                end
            end
        end
        %ieAddObject(scene{1}); ieAddObject(scene{2}); sceneWindow;
        
        % Compute optical images from the scene
        for ii = 1:2
            OIs{ii} = oiCompute(oi,scene{ii});
        end
        % ieAddObject(OIs{1}); ieAddObject(OIs{2}); oiWindow;
        
        ois = oiSequence(OIs{1}, OIs{2}, sampleTimes, modulation, ...
            'composition', composition);
        % ois.visualize;  % Not working right.  Something about image scale
        
        
    otherwise
        error('Unknown type %s\n',oisType);
end

%%