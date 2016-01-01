function scenes = rdtScenesLoad(varargin)
% Load multispectral scenes using Remote Data Toolbox
%
%   scenes = rdtScenesLoad(varargin)  
%
% Inputs:
%   varargin - name value pairs for the scene parameters
%  
%   'rdtConfigName' - repository name ('isetbio')
%   'fov'        - Field of view      (10)
%   'wave'       - Wavelength samples (400:10:680);
%   'nScenes'    - First n scenes     (1)
%   'sceneNames' - Specific scene names, these are the artifactId 
%    
% Outputs:
%   scenes - cell array of multispectral scenes
%
% Example:
%     scenes = rdtScenesLoad();
%     scenes = rdtScenesLoad('nScenes', 5); % Loads 5 scenes
%
% See also:
%   scarletScenesLoad
% 
% HJ/BW, VISTA TEAM, 2015

%% Parse input parameters
p = inputParser;
p.addParameter('nScenes', inf);
p.addParameter('wave', 400:10:680);
p.addParameter('rdtConfigName', 'isetbio');
p.addParameter('fov', 10);
p.parse();

nScenes   = p.Results.nScenes;
wave      = p.Results.wave;
rdtName   = p.Results.rdtConfigName;
fov       = p.Results.fov;

%% Init remote data toolbox client
rdt = RdtClient(rdtName);  % using rdt-config-scien.json configuration file
rdt.crp('/L3/faces');
files = rdt.listArtifacts();
nScenes = min(nScenes, length(files));

% load scene files
scenes = cell(nScenes, 1);
for ii = 1 : nScenes
    data = rdt.readArtifact(files(ii).artifactId);
    scene = sceneSet(data.scene, 'wave', wave); % adjust wavelength
    scenes{ii} = sceneSet(scene, 'h fov', fov);
end

end