%% Upload to RDT (optional)
%
% Make your cmosaic

% Then do this 
if rdtUploadFlag
    curDir = pwd;
    
    % We want to write a wrapper that puts, say, this script up into the
    % same artifact so we could reproduce the cone mosaic.
    chdir(isetbioRootPath,'local');
    save('rings-rays2.mat','cMosaic');
    filename1 = fullfile(isetbioRootPath,'local','rings-rays2.mat');
    client = RdtClient('isetbio');
    
    client.credentialsDialog();
    client.crp('/resources/data/cmosaics')
    version1 = '1';
    artifact = client.publishArtifact(filename1, 'version', version1);
    % client.openBrowser;
    
    chdir(curDir);
end

%% Download from RDT (optional)

if rdtDownloadFlag
    
    rdt = RdtClient('isetbio');
    rdt.crp('/resources/data/cmosaics');
    data = rdt.readArtifact('rings-rays', 'type', 'mat');
    cMosaic = data.cMosaic;

end

%% Save a cone mosaic data set related to vernier acuity (s_Layers)
% 
% The file vaConeMosaic.mat contains a cMosaic computed for a vernier
% acuity target with some eye movements.  Part of the WLVernierAcuity
% project.
%
% It also contains the computed cone photocurrent (alignedC) in response an
% aligned stimulus (with eye movements).

% Push it up to the site
rdt = RdtClient('isetbio');
rdt.credentialsDialog();

rdt.crp('/resources/data/cmosaics');
rdt.listArtifacts('type','mat');
filename1 = fullfile(isetbioRootPath,'local','vaConeMosaic.mat');
artifact = rdt.publishArtifact(filename1, 'version', '1');

%% To download the VA cone Mosaic data into the variable coneMosaicData

rdt = RdtClient('isetbio');
rdt.crp('/resources/data/cmosaics');
rdt.listArtifacts('type','mat','print',true);

data = rdt.readArtifact('vaConeMosaic', 'type', 'mat');
coneMosaicData = data.cMosaic;
alignedC       = data.alignedC;

%%

%%