%% s_coneAbsorbanceData
%
% These are the cone photopigment absorbances without any other
% ocular transmittance factors.

%% T_XX is in the Psychtoolbox-3
% Make sure that is on your path

foo = load('T_log10coneabsorbance_ss');

absorbance = foo.T_log10coneabsorbance_ss;
wave = SToWls(foo.S_log10coneabsorbance_ss);

%%
comment = 'Stockman/Sharpe cone absorbance from PTB.';
fname = fullfile(isetRootPath,'data','human','coneAbsorbance');

ieSaveSpectralFile(wave(:),absorbance',comment,fname);

%%

wave = 400:5:800;
a = ieReadSpectra('coneAbsorbance',wave);
vcNewGraphWin;
plot(wave,a)
hold on;

wave = 400:10:700;
a = ieReadSpectra('coneAbsorbance',wave);
plot(wave,a)
