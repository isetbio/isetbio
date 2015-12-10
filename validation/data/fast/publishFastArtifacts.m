%% One-time script used to move artifacts from local repo to Archiva.
%
% We have a new Archiva server hosted by AWS and the Brainard lab.  Browse
% or search the repositories with your web browser:
%   http://52.32.77.154/#browse
%   http://52.32.77.154/#search
%
% There is a repository called "isetbio".  It's located here:
%   http://52.32.77.154/#browse~isetbio
%   http://52.32.77.154/repository/isetbio/
%
% This script takes several validation artifacts currently part of the
% isetbio code repository and publishes them to our isetbio Archiva
% repository.
%
% Here is where I put things and how I chose to name them.  This all seemed
% reasonable an natural to me, but I'd be happy to revisit these choices.
%
% All artifacts go into the "isetbio" repository.
%
% All artifacts go into a parent group called "validation.fast".  There
% there is also a "validation.full", which was created from our Dropbox
% folder.
%
% Each folder in "isetbio/validation/data/fast" defines a subgroup.  So we
% have groups like "validation.fast.color", "validation.fast.wavefront",
% etc.
%
% Artifact names are taken from the original file names, minus the "v_"
% prefix and minus the "_FastGroundTruthDataHistory" suffix.  So for
% example, "v_Colorimetry_FastGroundTruthDataHistory.mat" gets the artifact
% name "Colorimetry"
%
% All artifacts have type "mat".
%
% Artifact versions come from the contents of each mat-file.  For example,
% "v_Colorimetry_FastGroundTruthDataHistory.mat" contains a single version
% of its data, called "run00001".  So "run00001" is the artifact version
% number.  If there were also a "run00002" this would be a published as
% another version of the same artifact.
%
% 27 November 2015
% benjamin.heasly@gmail.com

clear;
clc;

%% Find our data files.
thisScript = which('publishFastArtifacts.m');
validationFolder = fileparts(thisScript);

d = dir(validationFolder);
nDirs = numel(d);

dataFiles = [];

for ii = 1:nDirs
    if ~d(ii).isdir || '.' == d(ii).name(1)
        continue;
    end
    
    localFolder = fullfile(validationFolder, d(ii).name);
    remotePath = rdtFullPath({'', 'validation', 'fast', d(ii).name});
    
    matFiles = dir([localFolder '/*.mat']);
    nMatFiles = numel(matFiles);
    for jj = 1:nMatFiles
        % drop "_v" and "_FastGroundTruthDataHistory".
        nameParts = strsplit(matFiles(jj).name, '_');
        artifactId = nameParts{2};
        
        dataFile = struct( ...
            'subFolder', d(ii).name, ...
            'artifactId', artifactId, ...
            'localFolder', localFolder, ...
            'localFile', fullfile(localFolder, matFiles(jj).name), ...
            'remotePath', remotePath, ...
            'type', 'mat');
        
        dataFiles = [dataFiles, dataFile];
    end
end

%% Publish each version within each data file.
client = RdtClient('isetbio');
client.credentialsDialog();

nDataFiles = numel(dataFiles);
for ii = 1:nDataFiles
    dataFile = dataFiles(ii);
    data = load(dataFile.localFile);
    
    % expand version history into temp files
    dataVersions = fieldnames(data);
    for dv = dataVersions
        versionName = dv{1};
        dataVersion = data.(versionName);

        tempFolder = fullfile(tempdir(), 'isetbio-validation', dataFile.subFolder);
        if 7 ~= exist(tempFolder, 'dir')
            mkdir(tempFolder);
        end
        
        tempFile = fullfile(tempFolder, [dataFile.artifactId '-' versionName '.mat']);
        save(tempFile, '-struct', 'dataVersion');
        
        % pulish each temp file
        fprintf('Publishing %d of %d:\n', ii, nDataFiles)
        disp(tempFile)
        client.crp(dataFile.remotePath);
        client.publishArtifact(tempFile, ...
            'artifactId', dataFile.artifactId, ...
            'version', versionName);
    end
end
