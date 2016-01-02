%% s_rdtSummary
%
% Check on the SCIEN and ISETBIO Archiva repositories
%
% BW ISETBIO Team, Copyright 2016

ieInit

%% Summarize SCIEN

rdt = RdtClient('scien');
artifacts = rdt.listArtifacts;
rPaths = rdt.listRemotePaths;
fprintf('SCIEN repository contains\n %d artifacts in %d remote paths\n',length(artifacts),length(rPaths));

%% Summarize ISETBIO
rdt = RdtClient('isetbio');
artifacts = rdt.listArtifacts;
rPaths = rdt.listRemotePaths;
fprintf('ISETBIO repository contains\n %d artifacts in %d remote paths\n',length(artifacts),length(rPaths));

%%  Get one artifact from the isetbio repository

% There is a scene6 artifact that is a matlab file.  This is
% 'City view. Spatial and illuminant data available. Scene was recorded in
% the Picoto area, Braga, Minho reg?' 

% ** NOT WORKING **
try
    testA = rdt.searchArtifacts('scene6');
    data = rdt.readArtifact(testA.artifactId);
catch
    disp('Still not working');
end

% Error using gradleFetchArtifact (line 100)
% error status 1 (Starting a new Gradle Daemon for this build (subsequent builds will be faster).
% :fetchIt FAILED
% 
% FAILURE: Build failed with an exception.
% 
% * Where:
% Build file '/Users/wandell/Github/remoteDataToolbox/external/gradle/fetch.gradle' line: 28
% 
% * What went wrong:
% Execution failed for task ':fetchIt'.
% > Could not resolve all dependencies for configuration ':fetch'.
%    > Could not find any matches for :scene6:+ as no versions of :scene6 are available.
%      Searched in the following locations:
%          http://52.32.77.154/repository/isetbio//scene6/maven-metadata.xml
%          http://52.32.77.154/repository/isetbio//scene6/
%      Required by:
%          :gradle:unspecified
%          

%% Download via the URL works
tmp = tempname;
websave(tmp,testA.url);
load(tmp,'scene');
ieAddObject(scene); sceneWindow;
delete(tmp);

%%