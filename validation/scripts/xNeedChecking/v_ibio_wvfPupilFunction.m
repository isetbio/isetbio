%% v_wvfPupilFunction
%
% Calculate and explore the pupil function.
%
% This is the the complex function over the surface of the pupil that
% describes the wavefront aberration.  The function combines the
% aberrations from the measured Zernicke Polynomial coefficients and the
% amplitude imposed (later) by the Stiles Crawford Effect.
% 
% This script will explore each of the parts and their inter-relationship.
% Then we will hopefully convert the pupilfunction to the PSF in another,
% simple script.
%
% (c) Stanford VISTA Team, 2012 (bw)

%% 
s_initISET

%% Set values in millimeters
wvfParams0 = wvfCreate('measured pupil',6,'calculated pupil',3);

% Sample data
sDataFile = fullfile(wvfRootPath,'data','sampleZernikeCoeffs.txt');
theZernikeCoeffs = importdata(sDataFile);

nSubjects = size(theZernikeCoeffs,2);

% Possibly used for display.  Shouldn't be here, though.
nRows = ceil(sqrt(nSubjects));
nCols = ceil(nSubjects/nRows);

%% Show the diffraction limited case
theWavelength = wvfGet(wvfParams0,'wave');
wvfParams0 = wvfComputePupilFunction(wvfParams0);
vcNewGraphWin; imagesc(angle(wvfGet(wvfParams0,'pupil function',theWavelength)))

% Initialize Stiles Crawford
wvfParams0 = wvfSet(wvfParams0,'sce params',sceCreate(theWavelength,'none'));

% This is for 550 nm
thisSubject  = 1;
wvfParams = wvfSet(wvfParams0,'zcoeffs',theZernikeCoeffs(:,thisSubject));
wvfParams = wvfComputePupilFunction(wvfParams);
pupilF = wvfGet(wvfParams,'pupil function',theWavelength);
vcNewGraphWin;  mesh(angle(pupilF))
% Add the x,y units, which describe the pupil, I think, in some scale
% related to mm, but not sure which.

% Now for 650nm
theWavelength = 650;
wvfParams = wvfSet(wvfParams,'wave',theWavelength);
wvfParams = wvfSet(wvfParams,'sce params',sceCreate(theWavelength,'none'));
wvfParams = wvfComputePupilFunction(wvfParams);
pupilF = wvfGet(wvfParams,'pupil function',theWavelength);
vcNewGraphWin;  imagesc(angle(pupilF))

% Add a SCE
wvfParams = wvfSet(wvfParams,'sce params',sceCreate(theWavelength,'berendschot_data'));
wvfParams = wvfSet(wvfParams,'wave',theWavelength);
wvfParams  = wvfComputePupilFunction(wvfParams);
pupilF = wvfGet(wvfParams,'pupil function',theWavelength);

vcNewGraphWin;  imagesc(angle(pupilF))

% It looks like the center of the SCE effect is not in the center of pupil.
% Maybe that is right?  Don't really understand.
vcNewGraphWin;  mesh(abs(pupilF))


%% End