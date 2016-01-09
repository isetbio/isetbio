% v_rgcExternal
% 
% A validation for the coupled GLM model of retinal ganglion cells.
% This file compares finds the response of an ON Parsol RGC mosaic using
% both the isetbio implementation of the GLM and the Chichilnisky Lab's
% code. The ON Parasol RGC parameters were fit to physiological data from
% the Chichilnisky lab labeled as 2012-08-09-3.
% 
% The Remote Data Toolbox is used to download the necessary data from an
% Archiva server. Access is password protected; email the isetbio contact
% for access. Make sure you have the Remote Data Toolbox repository 
% (available on the isetbio github at isetbio/RemoteDataToolbox) added to
% your Matlab path. Four files are loaded using the RDT:
% 
%       1. testmovishort: The movie stimulus. In this case, a short 
%       sequence of natural scenes with simulated eye movements.
%       2. parasol_on_1205: The fitted parameters for the ON Parasol mosaic. 
%       This consists  of 39 separate "fittedGLM" structures.
%       ******** Right now, this is only 1 cell; will change soon.*******
%       3. pairspikesall: The measured spiking activity of the mosaic. This 
%       input is required for the coupled GLM.
%       4. xvalall: The output of the same simulation from code used by the
%       Chichilnisky Lab.
% 
% The script uses the osIdentity object for the cone outersegment
% representation, which simply copies the RGB data from the isetbio scene
% for each frame of the movie. This is due to the fact that the GLM
% is input-referred as opposed to cone current-referred. A new version of
% the isetbio GLM implementation will be cone current-referred and will be
% added in the near future.
% 
% (c) isetbio team
% JRG 12/2015
%% 
clear
ieInit;

%% Load stimulus movie

% % Load stimulus movie using RemoteDataToolbox
client = RdtClient('isetbio');
% client.credentialsDialog();
client.crp('resources/data/rgc')
[data, artifact] = client.readArtifact('testmovieshort', 'type', 'mat');
testmovieshort = data.testmovieshort;

% % Load movie from local data
% glmFitPath = '/Users/james/Documents/matlab/NSEM_data/';
% movieFiles = dir([glmFitPath 'testmovie*']);
% load([glmFitPath movieFiles(1).name], 'testmovie');
% testmovieshort = testmovie.matrix(:,:,1:1+5*120);
% % testmovieshort = (permute(testmovieshort, [2 1 3]));

%% Generate outersegment object
scene = 0; sensor = 0; % not needed because using osIdentity object

% osIdentity is used for input to the rgcPhys isetbio object
os2 = osCreate('identity');
os2 = osSet(os2, 'rgbData', double(testmovieshort));

%% Generate RGC object
% parameters not used because data is loaded from mat files
eyeAngle = 180; % degrees
eyeRadius = 10; % mm
eyeSide = 'right';

rgc2 = rgcPhys(scene, sensor, os2, eyeSide, eyeRadius, eyeAngle);
rgc2 = rgcSet(rgc2,'numberTrials',20);
rgc2 = rgcCompute(rgc2, os2);

rgc2linear = mosaicGet(rgc2.mosaic{1},'linearResponse');
rgc2psth = mosaicGet(rgc2.mosaic{1},'psthResponse');

%% Load validation data
% Load RDT version of output from the Chichilnisky Lab's code
client = RdtClient('isetbio');
% client.credentialsDialog();
client.crp('resources/data/rgc')
[data, artifact] = client.readArtifact('xvalall', 'type', 'mat');
xvalall = data.xvalall;

% % Load local copy
% load('xvalall.mat');

%% Compare isetbio output and Chichilnisky Lab output
vcNewGraphWin;
for i = 1%:36
    % Measure difference between outputs
    minlen = min([length(rgc2psth{i}) length(xvalall{i}.psth)]);
    diffpsth(i) = sum(abs(rgc2psth{i}(1:minlen) - xvalall{i}.psth(1:minlen)))./sum(.5*(rgc2psth{i}(1:minlen) + xvalall{i}.psth(1:minlen)));
    
    % % Plot difference
    % subplot(6,7,i); hold on;
    % plot(rgc2psth{i}(1:minlen)-xvalall{i}.psth(1:minlen),'b','linewidth',1);

    % Plot output of isetbio code
    plot(rgc2psth{i}(1:minlen),'r ','linewidth',3);
    hold on;
    % Plot output of Chichilnisky Lab code
    plot(xvalall{i}.psth(1:minlen),':k','linewidth',2);

    axis([0 6285 0  100]);%max(rgc2psth{i}(1:minlen))])
    [maxv, maxi] = max(rgc2psth{i}(1:minlen)-xvalall{i}.psth(1:minlen)); 
    title(sprintf('maxv = %.1f, maxi = %d',maxv,maxi));
end

% vcNewGraphWin; plot(diffpsth,'x')

%% Compare linear outputs
% figure;
% for i = 1:39
%     subplot(6,7,i);
% %     plot(rgc2psth{i});
% %     hold on;
% %     plot(xvalall{i}.psth);
%     
%     plot(rgc2linear{i}(1:minlen/10)-xvalall{i}.lcif_const(1:10:minlen),'b','linewidth',2);
%     hold on;
% %     plot(xvalall{i}.cif0);
%     
% end
