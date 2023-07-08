%% v_ibio_wvfStilesCrawford

%% Some informative text
% UnitTest.validationRecord('SIMPLE_MESSAGE', 'Validate wavefront Stiles-Crawford Effect code.');

%%
ieInit;

%% For plotting limits
maxMIN = 2; maxMM  = 1; waveIdx = 1;
thisWave = 550;

% Set up wvf structure and sce structure
wvf = wvfCreate;
wvf = wvfSet(wvf,'zcoeffs',[0.2 0.75],{'defocus', 'oblique_astigmatism'});
sceP = sceCreate(thisWave,'berendschot_data','centered');

%% No Stiles Crawford effect

wvf = wvfSet(wvf,'sce params',[]);
wvf = wvfCompute(wvf);
sce1DFig2 = ieNewGraphWin; hold on
wvfPlot(wvf,'1d psf angle','min',[],maxMIN,'no window');
hold on; grid on; 

%% Include the SCE 

wvfSCE = wvfSet(wvf,'sce params',sceP);
wvfSCE = wvfCompute(wvfSCE,'compute sce',true);
[f,p] = wvfPlot(wvfSCE,'1d psf angle','min',[],maxMIN,'no window');
set(p,'color','b','LineStyle','--')
legend({'No SCE','Yes SCE'});

%% Show the two aperture functions

aperture = wvfGet(wvf,'aperture');
apertureSCE = wvfGet(wvfSCE,'aperture');
pupilPos = wvfGet(wvf,'pupil positions',thisWave,'mm');
ieNewGraphWin([],'wide');
tiledlayout(1,2);
nexttile; mesh(pupilPos,pupilPos,aperture); title('No SCE');
xlabel('Position (mm)'); ylabel('Position (mm)');
nexttile; mesh(pupilPos,pupilPos,apertureSCE); title('SCE');

%% END


