%% Data for Pokorny Smith lens density as a function of age
%
% Aging of the human lens
% Table I:   Pokorny Smith and Lutze, 1987
%
% Tested in Optical density of the human lens, Xu, Pokorny, Smith, 1997
%

%% Wavelength

w = 400:10:650;

% Two components
TL1 = [
    0.6000
    0.5100
    0.4330
    0.3770
    0.3270
    0.2950
    0.2670
    0.2330
    0.2070
    0.1870
    0.1670
    0.1470
    0.1330
    0.1200
    0.1070
    0.0930
    0.0800
    0.0670
    0.0530
    0.0400
    0.0330
    0.0270
    0.0200
    0.0130
    0.0070
    0];

TL2 = [
    1.0000
    0.5830
    0.3000
    0.1160
    0.0330
    0.0050
    0
    0
    0
    0
    0
    0
    0
    0
    0
    0
    0
    0
    0
    0
    0
    0
    0
    0
    0
    0];

TYoung = [
    1.7649
    1.1374
    0.7240
    0.4876
    0.3413
    0.2629
    0.2279
    0.2046
    0.1834
    0.1675
    0.1537
    0.1378
    0.1230
    0.1102
    0.0986
    0.0859
    0.0742
    0.0615
    0.0488
    0.0381
    0.0297
    0.0223
    0.0170
    0.0117
    0.0053
    0.0032];

%% Check
ieNewGraphWin;
plot(w,TL1, 'o',w,TL2,'x');

data = [TL1(:), TL2(:), TYoung(:)];
fname = fullfile(isetbioRootPath,'data','datafiles','human','lens','lensPokornySmith.mat');
fullpathname = ieSaveSpectralFile(w, data, 'Model from Pokorny, Smith, Lutze, Table I',fname);

%% Check

extraWave = 400:10:700;
foo = ieReadSpectra(fullpathname,(400:10:700));
plot(extraWave,foo);
hold on;
plot(extraWave,foo(:,1) + foo(:,2));

%% END
