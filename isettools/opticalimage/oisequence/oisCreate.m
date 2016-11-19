function [ois, varargout] = oisCreate(oisType,varargin)
% Create an oiSequence
%
%     ois = oisCreate(oisType,param/val)
%
% Examples
%
%  Harmonic -  TODO:  We can add parameters to the params(:) struct that
%  specify the scene sets.
%   
%    params(1) = harmonicP; params(2) = params(1); params(2).contrast = 0;
%    [ois, scenes] = oisCreate('harmonic','params',params);
%    ois.visualize;
%
%  Vernier - TODO
%
% BW ISETBIO Team, 2016

%%
p = inputParser;
p.addRequired('oisType',@ischar)
p.addParameter('params',[],@isstruct);
p.parse(oisType,varargin{:});

oisType = ieParamFormat(p.Results.oisType);

% Many of the types have a params structure that we will pass along
params = p.Results.params;

%%

switch oisType
    case 'harmonic'
        % oisCreate('harmonic',params);
        % params:  A two element array as returned by harmonicP;
        if length(params) ~= 2, error('Two param sets needed.'); end
        scene = cell(1,2);
        OIs = cell(1, 2);
        
        % We need a way to send in scene parameters for setting
        imgFov = .5 ;      % image field of view
        vDist  = 0.3;          % viewing distance (meter)
        for ii=1:2
            scene{ii} = sceneCreate('harmonic',params(ii));
            scene{ii} = sceneSet(scene{ii}, 'h fov', imgFov);
            scene{ii} = sceneSet(scene{ii}, 'distance', vDist);
            scene{ii} = sceneSet(scene{ii},'name',...
                sprintf('F %d C %.2f', params(ii).freq, params(ii).contrast));
        end
        % ieAddObject(scene{1}); ieAddObject(scene{2}); sceneWindow;
                
        % Compute optical images from the scene
        oi = oiCreate('wvf human');
        for ii = 1 :2
            OIs{ii} = oiCompute(oi,scene{ii});
        end
        
        %% Build the oiSequence
        
        % We need a way to send in the weights (modulationFunction)
        
        % We build the stimulus using a time series of weights. We have the mean
        % field on for a while, then rise/fall, then mean field.
        zTime = 50;   % Mean field beginning and end (ms)
        stimWeights = fspecial('gaussian',[1,50],15);
        stimWeights = ieScale(stimWeights,0,1);
        weights = [zeros(1, zTime), stimWeights, zeros(1, zTime)];
        
        % Temporal samples.  Typically 1 ms, which is set by the parameter in the
        % cone mosasic integration time.  That time is locked to the eye movements.
        tSamples = length(weights);
        sampleTimes = 0.002*(1:tSamples);  % Time in sec
        % vcNewGraphWin; plot(1:tSamples, weights,'o');
        % xlabel('Time (ms)');
        
        % The weights define some amount of the constant background and some amount
        % of the line on the same constant background
        ois = oiSequence(OIs{2}, OIs{1}, ...
            sampleTimes, weights, ...
            'composition', 'blend');
        
        varargout{1} = scene;
        
    otherwise
        error('Unknown type %s\n',oisType);
end

%%