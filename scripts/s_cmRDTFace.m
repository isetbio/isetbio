%% Make a cone mosaic data of faces and put it up on RDT
%
% s_cmRDTFace
%
% See rdt.removeArtifact proposal below
%
% BW ISETBIO Team, 2017

%% Make a cone mosaic of a face image

scene = sceneFromFile(fullfile(isetRootPath,'data','images','rgb','faceMale.jpg'),'rgb');
scene = sceneSet(scene,'fov',1.2);
scene = sceneAdjustIlluminant(scene,'D65.mat');
vcAddObject(scene); sceneWindow;

oi = oiCreate;
oi = oiCompute(oi,scene);
vcAddObject(oi); oiWindow;

%%
cMosaic = coneMosaic;
cMosaic.setSizeToFOV(sceneGet(scene,'fov'));
cMosaic.emGenSequence(60);

disp('Computing cone mosaic current');
cMosaic.compute(oi);
cMosaic.computeCurrent(oi);
cMosaic.window;

%%
faceFile = fullfile(isetbioRootPath,'local','coneMosaicDataFace.mat');
save(faceFile,'cMosaic');

rdt = RdtClient('isetbio');
rdt.credentialsDialog();  % wandell, Jxxx4XX

% You can change other places.
rdt.crp('/resources/data/cmosaics')
rdt.listArtifacts('recurseive',true,'print',true);
version1 = '1';
rdt.publishArtifact(faceFile, 'version', version1);
a = rdt.listArtifacts('print',true);


%%  If you want to update the face mosaic try this

% rdt = RdtClient('isetbio');
% rdt.credentialsDialog();  % wandell, Jxxx4XX
% rdt.crp('/resources/data/cmosaics')
% 
% a = rdt.listArtifacts('print',true);
% 
% % Choose the face mosaics
% % facemosaics = a(1:2);
% 
% % BSH: let's look specifically at the image-artifact version 1 folder
% % it should contain the jpg and various metadata files
% % artifactFolderHack = rdtArtifact('url', fileparts(a.url));
% % rd.openBrowser(artifactFolderHack);
% 
% % Remove the file
% 
% % We should add a method to RdtClient called 'removeArtifacts'
% %  It should do this
% deleted = rdtDeleteArtifacts(rdt.configuration, facemosaics, 'allFiles', true);
% 
% % This returns an empty artifact, which is good
% a = rd.listArtifacts('print',true);

%%

