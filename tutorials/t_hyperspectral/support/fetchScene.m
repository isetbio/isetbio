function theScene = fetchScene(sceneNo, varargin)
% Load a hyperspectral scene, either from the remote server or from the disk
%
% Syntax:
%    theScene = fetchScene(sceneNo, varargin)
%
% Description:
%    Load a hyperspectral scene, stored either on the Brainard Lab archiva
%    server (remote), or on the disk (locally)
%
% Inputs:
%    sceneNo     - the scene no, 1, 2, ...
%    database    - the hyperspectral data base identifier
%    destination - where to save the downloaded scene
%
% Outputs:
%    theScene    - an ISETBio scene generated from the retrieved 
%                  hyperspectral image data
%
% Optional key/value pairs:
%    origin      - a string set to either 'remote' (for remote scene
%                 fetching) or to some path to the scene, .e.g.,
%                 '~/Desktop'. In that case the scene filename must be
%                 'sceneN.mat', where N is the passed sceneNo
%
% 7/24/18  npc  Wrote it
%
%{
    % Load scene 3 from the remote server
    sceneNo = 3;
    theScene = fetchScene(sceneNo, ...
        'origin', 'remote',  ...
        'database', 'manchester_database/2002');
%}

%{
    % Load scene 3 from the remote server
    sceneNo = 3;
    theScene = fetchScene(sceneNo, ...
        'origin', 'remote',  ...
        'database', 'manchester_database/2004');
%}

%{
    % Load scene 3 from ~/Desktop/images/scene3.mat
    sceneNo = 3;
    theScene = fetchScene(sceneNo, 'origin', '~/Desktop/images');
%}


    p = inputParser;
    p.addParameter('origin', 'remote');
    p.addParameter('database', 'manchester_database/2004');
    p.addParameter('destination', './resources');
    p.parse(varargin{:});
    origin = p.Results.origin;
    database = p.Results.database;
    destination = p.Results.destination;
    
    if (strcmp(origin, 'remote'))
        % Set up remote data toolbox client
        rd = RdtClient('isetbio');
        % Look in the Manchester 2004 data set
        rd.crp(fullfile('/resources/scenes/hyperspectral', database));
        % Get a list of all the artifact names
        a = rd.listArtifacts('type','mat');
        if (numel(a) < sceneNo)
            error('There are only %d scenes in this data set. You specified scene #%d', numel(a), sceneNo)
        end
        % Load the hyperspectral scene
        sData = rd.readArtifact(a(sceneNo),'type','mat');
        % Generate an isetbio scene from the hyperspectral data
        scene = sceneFromBasis(sData);
        % Save locally
        saveFileName = fullfile(destination, database, sprintf('scene%d', sceneNo));
        save(saveFileName, 'scene');
        fprintf('Scene saved in %s\n', saveFileName);
    else
        localSceneFileName = fullfile(origin, database, sprintf('scene%d', sceneNo));
        load(localSceneFileName,'scene');
    end
    theScene = scene;
end