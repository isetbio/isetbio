function v_ibio_sceneHCCompress(varargin)
%
% Validate hyperspectral scene data and compression
%
% Uses linear model compression.
%
% Copyright ImagEval Consultants, LLC, 2012
% 
%     varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
% end
% 
% 
% %% Here is the actual code
% function ValidationFunction(runTimeParams)

%% Read in the scene
fName = fullfile(isetRootPath,'data','images','multispectral','StuffedAnimals_tungsten-hdrs');
scene = sceneFromFile(fName,'multispectral');
UnitTest.validationData('sceneA', scene);
    
% Have a look at the image
% Show scene in window
sceneWindow(scene);

% Plot the illuminant
scenePlot(scene,'illuminant photons');

%% Compress the hypercube requiring only 95% of the var explained
vExplained = 0.95;
[imgMean, imgBasis, coef] = hcBasis(sceneGet(scene,'photons'),vExplained);

%% Save the data 
wave        = sceneGet(scene,'wave');
basis.basis = imgBasis;
basis.wave  = wave;

comment = 'Compressed using hcBasis with imgMean)';

illuminant = sceneGet(scene,'illuminant');
% illuminant.wavelength = scene.spectrum.wave;
% illuminant.data = scene.illuminant.data;
oFile = fullfile(isetbioRootPath,'local','deleteMe.mat');
ieSaveMultiSpectralImage(oFile,coef,basis,comment,imgMean,illuminant);

%% read in the data
wList = 400:5:700;
scene2 = sceneFromFile(oFile ,'multispectral',[],[],wList);

% It is very desaturated
sceneWindow(scene2);


%% Now require that most of the variance be plained
vExplained = 0.99;
[imgMean, imgBasis, coef] = hcBasis(sceneGet(scene,'photons'),vExplained);
fprintf('Number of basis functions %.0f\n',size(imgBasis,2));

%% Save the data 
wave        = sceneGet(scene,'wave');
basis.basis = imgBasis;
basis.wave  = wave;

comment = 'Compressed using hcBasis with imgMean)';

illuminant = sceneGet(scene,'illuminant');
% illuminant.wavelength = scene.spectrum.wave;
% illuminant.data = scene.illuminant.data;
ieSaveMultiSpectralImage(oFile,coef,basis,comment,imgMean,illuminant);

%% read in the data
wList = 400:5:700;
scene2 = sceneFromFile(oFile ,'multispectral',[],[],wList);
sceneWindow(scene2);

%% Clean up the temporary file.
delete(oFile);

end
