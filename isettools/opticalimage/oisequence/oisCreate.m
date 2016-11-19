function [ois, varargout] = oisCreate(oisType,composition, modulation, varargin)
% Create an oiSequence
%
%     ois = oisCreate(oisType,param/val)
%
% Examples
%
%  Harmonic -  TODO:  We can add parameters to the params(:) struct that
%  specify the scene sets.
%   
% We need a way to send in the weights (modulationFunction)
%
% We build the stimulus using a time series of weights. We have the mean
% field on for a while, then rise/fall, then mean field.
%
% stimWeights = ieScale(fspecial('gaussian',[1,50],15),0,1);
% weights = [zeros(1, 30), stimWeights, zeros(1, 30)];
% 
% Example:
%
%    hparams(1) = harmonicP; 
%    hparams(2) = hparams(1); hparams(2).contrast = 0;
%    sparams.fov = 0.3; 
%    [ois, scenes] = oisCreate('harmonic','blend',weights, 'hparams',hparams,'sparams',sparams);
%    ois.visualize;
%
%  Vernier - TODO
%
% BW ISETBIO Team, 2016

%% Inputs
p = inputParser;
p.addRequired('oisType',@ischar)
p.addRequired('composition',@ischar);
p.addRequired('modulation');

p.addParameter('sampleTimes',[],@isvector);
p.addParameter('hparams',[],@isstruct);
p.addParameter('sparams',[],@isstruct);

p.parse(oisType,composition,modulation,varargin{:});

oisType = ieParamFormat(p.Results.oisType);
composition = p.Results.composition;
modulation  = p.Results.modulation;

sampleTimes = p.Results.sampleTimes;
if isempty(sampleTimes)
    sampleTimes = 0.001*((1:length(modulation))-1);
end

% Many of the types have a params structure that we will pass along
hparams = p.Results.hparams;   % Harmonics parameters (two structs needed)
sparams = p.Results.sparams;   % Scene parameters, only one set.

%%
switch oisType
    case 'harmonic'
        % oisCreate('harmonic',params);
        % params:  A two element array as returned by harmonicP;
        if length(hparams) ~= 2, error('Specify two harmonic param sets.'); end
        scene = cell(1,2);
        OIs = cell(1, 2);

        % Create basic harmonics
        for ii=1:2
            scene{ii} = sceneCreate('harmonic',hparams(ii));
            sname = sprintf('F %d C %.2f', hparams(ii).freq, hparams(ii).contrast);
            scene{ii} = sceneSet(scene{ii},'name',sname);
        end
        % ieAddObject(scene{1}); ieAddObject(scene{2}); sceneWindow;
        
        % Adjust both scenes based on sparams.
        fields = fieldnames(sparams);
        if ~isempty(fields)
            for ii=1:length(fields)
                for jj=1:2
                    val = eval(['sparams.',fields{ii}]);
                    scene{jj} = sceneSet(scene{jj}, fields{ii},val);
                end
            end
        end
        
        % Compute optical images from the scene
        oi = oiCreate('wvf human');
        for ii = 1:2
            OIs{ii} = oiCompute(oi,scene{ii});
        end
        % ieAddObject(OIs{1}); ieAddObject(OIs{2}); oiWindow;

        %% Build the oiSequence
        
        % The weights define some amount of the constant background and some amount
        % of the line on the same constant background
        ois = oiSequence(OIs{2}, OIs{1}, sampleTimes, modulation, ...
            'composition', composition);
        
        % Return the cell array of scenes.
        varargout{1} = scene;
        
    otherwise
        error('Unknown type %s\n',oisType);
end

%%