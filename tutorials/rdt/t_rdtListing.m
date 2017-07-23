%% List some of the remote data
%
% We list some of the paths (directories)
% We also list some of the artifacts (files)
%
% This script should continue running, although I think we need to move the
% data around to organize the stimuli, cone mosaics, and whatever else is
% there.
%
% BW, ISETBIO Team, 2017

%%
ieInit;

%% Get data from the RDT site
rdt = RdtClient('isetbio');

%% List only one directory deep 
rdt.crp('/resources');
rdt.listRemotePaths('print',true);  % 'all' false limits

%% Here are the hyperspectral scene directories
rdt.crp('/resources/scenes/hyperspectral');
rdt.listRemotePaths('print',true);

%% 
rdt.crp('/resources/data');
rdt.listRemotePaths('print',true);

%% These are coneMosaics stored for computations

rdt.crp('/resources/data/cmosaics');
rdt.listArtifacts('print',true);

%% These are stimuli stored for computations

rdt.crp('/resources/data/istim');
rdt.listArtifacts('print',true);

%%