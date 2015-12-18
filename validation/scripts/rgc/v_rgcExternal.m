clear

glmFitPath = '/Users/james/Documents/matlab/NSEM_data/';
movieFiles = dir([glmFitPath 'testmovie*']);
load([glmFitPath movieFiles(1).name], 'testmovie');
testmovieshort = testmovie.matrix(:,:,1:121+5*120);
% testmovieshort = (permute(testmovieshort, [2 1 3]));

scene = 0; sensor = 0; % not needed because using osIdentity object

os2 = osCreate('identity');
os2 = osSet(os2, 'rgbData', double(testmovieshort));

eyeAngle = 180; % degrees
eyeRadius = 3; % mm
eyeSide = 'right';

rgc2 = rgcPhys(scene, sensor, os2, eyeSide, eyeRadius, eyeAngle);
rgc2 = rgcSet(rgc2,'numberTrials',20);
rgc2 = rgcCompute(rgc2, os2);

