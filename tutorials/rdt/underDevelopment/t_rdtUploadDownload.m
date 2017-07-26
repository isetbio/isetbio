%% Upload to RDT (optional)
%
% We are using the Rdtrdt software on an AWS site to manage shared data.
% BW hopes that in the future we will be able to use Flywheel and the
% scitran rdt in order to simplify his life.
%
% But we currently have a lot of things up on the archiva site, and it will
% be with us for at least another year.  

% So, here are some examples of putting data on the site, listing the data
% in a directory, and getting data down from the site.
%
% See the ISETBIO wiki and the remote data toolbox wikie for more
% information. 
%   
%   https://github.com/isetbio/isetbio/wiki/ISETBIO-Data
%   https://github.com/isetbio/RemoteDataToolbox/wiki
%
% BW, ISETBIO Team, 2017.

% Then do this 
if rdtUploadFlag
    curDir = pwd;
    
    % We want to write a wrapper that puts, say, this script up into the
    % same artifact so we could reproduce the cone mosaic.
    chdir(isetbioRootPath,'local');
    % save('rings-rays2.mat','cMosaic');
    
    % This should be a file that have stored in the local directory.
    filename1 = fullfile(isetbioRootPath,'local','rings-rays2.mat');
    
    rdt = Rdtrdt('isetbio');
    rdt.credentialsDialog();
    
    % You can change other places.
    rdt.crp('/resources/data/cmosaics')
    version1 = '1';
    artifact = rdt.publishArtifact(filename1, 'version', version1);

    % To remove the artifact, use 
    
    
    % rdt.openBrowser;
    
    chdir(curDir);
end

%% Download from RDT (optional)

if rdtDownloadFlag
    
    rdt = Rdtrdt('isetbio');
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
rdt = Rdtrdt('isetbio');
rdt.credentialsDialog();

rdt.crp('/resources/data/cmosaics');
rdt.listArtifacts('type','mat');
filename1 = fullfile(isetbioRootPath,'local','vaConeMosaic.mat');
artifact = rdt.publishArtifact(filename1, 'version', '1');

%% To download the VA cone Mosaic data into the variable coneMosaicData

rdt = Rdtrdt('isetbio');
rdt.crp('/resources/data/cmosaics');
rdt.listArtifacts('type','mat','print',true);

data = rdt.readArtifact('vaConeMosaic', 'type', 'mat');
coneMosaicData = data.cMosaic;
alignedC       = data.alignedC;

%%

%%