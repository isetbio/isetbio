% v_wvfSVNVer121TestData
%
% Tests that we can reconstruct PSFs computed using version 121 of the
% toolbox (BrainardLab SVN server).  This was pretty early in the
% development, April 30, 2011, and prior to various major organizational
% changes.
%
% The version that generated the test data is saved as a branch on the
% BrainardLab SVN server, with the new test program that writes out the data. 
% If we need more test cases, can generate them there.
%
% This is not really an independent test of the underlying code, but does
% serve to verify that the basic computations still work as they did
% when we first got the code from Heidi Hofer.
%
% As of 7/4/12, it appears that there has been a switch in the spatial
% coordinate system of the PSFs. 
% I believe this happened because Heidi's code's x/y convention was mismatched
% to Matlab's row/col convention, and one of Brian's students switched this.  In
% addition, the y coordinate in Heidi's code ran in the opposite direction
% from the y coordinate in our code.  Again, we think current
% code matches the OSA conventions, based on a first principles reading of 
% the code as well as our current agreement with the PSF pictures in Autrusseau
% et al.
%
% As of 7/29/12, it also appears that LCA was getting added in with the wrong
% sign in the ver 121 code. At least, it is getting added in with the opposite
% sign in that code.
%
% As of 7/29/12, we have changed to pass the j=0 coefficient.  So we need to
% handle that here too.
%
% 7/4/12  dhb  Wrote first draft.
%
% (c) Wavefront Toolbox Team, 2012

%% Initialize
s = which('v_wvfSVNVer121TestData');
cd(fileparts(s));
clear; close all;
%s_initISET;

%% Plot diffraction (tends to compress the scale of observer calcs)
PLOT_DIFFRACTION = 0;

%% Check PSFs off the middle row/col
rowOffset = 5;

%% Load in a test data set computed with SVN Version 121 (BrainardLab server)
%
% This is a very early version of the code as we got it from Heidi Hofer
% and serves as a test that we haven't munged anything up since that point.
theFiles = {...
    'SVNVer121_subj1_calcpupil3_defocus0_wavelength550_sce0', ...
    'SVNVer121_subj4_calcpupil3_defocus0_wavelength550_sce0', ...
    'SVNVer121_subj1_calcpupil3_defocus0_wavelength450_sce0', ...
    'SVNVer121_subj1_calcpupil5_defocus0_wavelength550_sce0', ...
    'SVNVer121_subj1_calcpupil3_defocus2_wavelength550_sce0', ...
    'SVNVer121_subj1_calcpupil3_defocus0_wavelength550_sce1', ...
    'SVNVer121_subj1_calcpupil3_defocus1_wavelength550_sce1', ...
    };

for i = 1:length(theFiles)
    testData = load(fullfile(wvfRootPath,'data','ver121data',theFiles{i}));
    
    %% Set up parameters structure to match old code, for diffraction limited computation
    wvf0 = wvfCreate;
    wvf0 = wvfSet(wvf0,'measured pupil size',testData.measpupilMM);
    wvf0 = wvfSet(wvf0,'measured wl',testData.nominalFocusWavelength);
    wvf0 = wvfSet(wvf0,'spatial samples',testData.sizeOfFieldPixels);
    wvf0 = wvfSet(wvf0,'ref pupil plane size',testData.sizeOfFieldMM);
    wvf0 = wvfSet(wvf0,'calc pupil size',testData.calcpupilMM);
    wvf0 = wvfSet(wvf0,'calc wavelengths',testData.theWavelength);
    if (testData.DOSCE == 1)
        sce = sceCreate(testData.theWavelength,'berendschot_data');
        wvf0 = wvfSet(wvf0,'sce params',sce);
    else
        sce = sceCreate(testData.theWavelength,'none');
        wvf0 = wvfSet(wvf0,'sce params',sce);
    end
    
    % Flip sign of defocus correction, because it was different in
    % version 121.
    lcaDiopters = wvfLCAFromWavelengthDifference(testData.nominalFocusWavelength,testData.theWavelength);
    wvf0 = wvfSet(wvf0,'calc observer focus correction',2*lcaDiopters+testData.defocusDiopters);
    
    % Compute diffraction limited PSF our way
    wvf0 = wvfSet(wvf0,'zcoeffs',zeros(size(testData.theZcoeffs)));
    wvf0 = wvfComputePSF(wvf0);
    diffracpsf0 = wvfGet(wvf0,'psf',testData.theWavelength);
    diffracpsfLine0 =  diffracpsf0(wvfGet(wvf0,'middle row')+rowOffset,:);
    diffracpsfLine0Centered = wvfGet(wvf0,'1d psf',testData.theWavelength);
    diffracarcmin0 = wvfGet(wvf0,'psf angular samples','min',testData.theWavelength);
    arcminperpix0 = wvfGet(wvf0,'psf arcmin per sample',testData.theWavelength);
    
    % Compute observer PSF our way.  We have changed so that we pass the
    % j = 0 coefficient, so need to prepend that here.
    wvf0 = wvfSet(wvf0,'zcoeffs',[0 ; testData.theZcoeffs]);
    wvf0 = wvfComputePSF(wvf0);
    psf0 = wvfGet(wvf0,'psf',testData.theWavelength);
    psfLine0 = psf0(wvfGet(wvf0,'middle row')+rowOffset,:);
    psfLineCentered0 = wvfGet(wvf0,'1d psf',testData.theWavelength);
    arcmin0 = wvfGet(wvf0,'psf angular samples','min',testData.theWavelength);
    
    % Make a comparison plot
    %
    % Notice that I need to pull out the column from the old
    % computation to match the row of the current computation, and
    % flip the sign of the row offset (which becomes a column offset
    % in the old code.)
    % I think the y direction is also inverted in the old calculation
    % but this turns out not to matter because the middle row corresponds
    % to y = 0.
    vcNewGraphWin; hold on
    if (PLOT_DIFFRACTION)
        plot(diffracarcmin0,diffracpsfLine0,'r','LineWidth',3);
        plot(testData.arcminutes,testData.diffracPSF(:,testData.whichRow-rowOffset),'g');
    end
    plot(arcmin0,psfLine0,'b','LineWidth',3);
    plot(testData.arcminutes,testData.thePSF(:,testData.whichRow-rowOffset),'g');
    plot(arcmin0,psfLineCentered0,'k','LineWidth',3);
    plot(testData.arcminutes,testData.thePSFCentered(:,testData.whichRow),'g');
    xlabel('arc minutes');
    ylabel('psf');
    title(LiteralUnderscore(theFiles{i}));
end

