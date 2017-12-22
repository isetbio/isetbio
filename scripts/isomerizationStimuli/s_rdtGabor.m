% Make a collection of Gabor isomerization stimuli
%
% Save the stimuli on the RDT for later use
%
% I would like to have this script be self-documenting so that we could
% have a web-page that described the stimuli on the RDT and how to download
% them.

% 

%%
ieInit

%%
chdir(tempdir);
if ~exist('Gabor','dir')
    mkdir('Gabor')
end
chdir('Gabor')
sDir = pwd;

%% Set up the frequency, contrast, and orientation directions

fList = logspace(-.5,1.6,8);  %8
fList = round(fList*100)/100; 
cList = logspace(-2,-0.3,5);  %5
cList = round(cList*100)/100;
aList = [0, pi/4, pi/2, 3*pi/4];
aListDeg = ieRad2deg(aList);

%% Create and save the iStim versionn of the Gaboor patch
pG.nSteps = 100; 
pG.GaborFlag = 0.2; 
pG.fov = 1;
pG.row = 160;
pG.col = 160;

% ff = 4; cc = 5; aa = 2;
for ff=1:length(fList)
    for cc = 1:length(cList)
        for aa = 1:length(aList)            
            fname = sprintf('Gabor-freq-%.2f-contrast-%.2f-orient-%d',fList(ff),cList(cc),aListDeg(aa));
            pG.contrast = cList(cc);
            pG.ang = aList(aa);
            pG.freq = fList(ff);
            iStim = ieStimulusGabor(pG);
            
            % Save the iStim in the save directory with the right name
            save(fullfile(sDir,[fname,'.mat']),'iStim'); 

            % vcAddObject(iStim.scene); sceneWindow;
            % vcAddObject(iStim.absorptions); sensorWindow;
            % coneImageActivity(iStim.absorptions,'dFlag',true);
            % sceneShowImage(iStim.scene);
        end
    end
end

% load('Gabor-freq-1.26-contrast-0.50-orient-90.mat');
% coneImageActivity(iStim.absorptions,'dFlag',true);

%% Upload to the RDT


%% A few for testing the Rieke biophysical data 
% These need to be temporally sampled by the sensor at half millisecond.
% So, we set pG.

%%
