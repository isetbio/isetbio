% CompareMonCalOverTime
%
% Compare two calibrations of a display.
%
% Refits all files to with common number of primary bases and gamma fitting method.
% These are prompted for.
%
% 1/20/05	dhb, bx		Wrote it.
% 2/12/10   dhb         Don't ask for load code, just prompt for name.
%           dhb         Better plots.  And, ask for which times to compare.
% 2/15/10   dhb         Fix input, not a string.
% 3/1/10    dhb         Allow drawing from different files, refitting data, etc.
% 6/23/11   dhb         Make a chromaticity plot of the comparison as well.
% 7/03/14   npc         Modifications for accessing calibration data using @CalStruct objects.
%  
%% Clear and close
clear; close all;

%% Get first calibration file and extract desired calibration
defaultFileName = 'LCDScreen';
thePrompt = sprintf('Enter first calibration filename [%s]: ',defaultFileName);
thenFileName = input(thePrompt,'s');
if (isempty(thenFileName))
    thenFileName = defaultFileName;
end
fprintf(1,'\nLoading from %s.mat\n',thenFileName);
[cal,cals] = LoadCalFile(thenFileName);
fprintf('Calibration file %s read\n',thenFileName);

% Print out available dates
fprintf('Calibration file contains %d calibrations\n',length(cals));
fprintf('Dates:\n');
for i = 1:length(cals)
    fprintf('\tCalibration %d, date %s\n',i,cals{i}.describe.date);
end

% Get which to compare
defaultThen = length(cals)-1;
thenIndex = input(sprintf('Enter number of first calibration to compare [%d]: ',defaultThen));
if (isempty(thenIndex))
    thenIndex = defaultThen;
end
if (thenIndex < 1 || thenIndex > length(cals))
    error('Calibration number out of range\n');
end
calThen = cals{thenIndex};

%% Get second calibration file and extract desired calibration.
% This can be the same file, or a different one.
defaultFileName = thenFileName;
thePrompt = sprintf('\nEnter second calibration filename [%s]: ',defaultFileName);
nowFileName = input(thePrompt,'s');
if (isempty(nowFileName))
    nowFileName = defaultFileName;
end
fprintf(1,'\nLoading from %s.mat\n',nowFileName);
[cal,cals] = LoadCalFile(nowFileName);
fprintf('Calibration file %s read\n',nowFileName);

% Print out available dates
fprintf('Calibration file contains %d calibrations\n',length(cals));
fprintf('Dates:\n');
for i = 1:length(cals)
    fprintf('\tCalibration %d, date %s\n',i,cals{i}.describe.date);
end

defaultNow = length(cals);
nowIndex = input(sprintf('Enter number of second calibration to compare [%d]: ',defaultNow));
if (isempty(nowIndex))
    nowIndex = defaultNow;
end
if (nowIndex < 1 || nowIndex > length(cals))
    error('Calibration number out of range\n');
end
calNow = cals{nowIndex};


% Specify @CalStruct objects that will handle all access to the calibration data.
[calThenOBJ, inputArgIsACalStructOBJ] = ObjectToHandleCalOrCalStruct(calThen);
clear 'calThen';

[calNowOBJ, inputArgIsACalStructOBJ] = ObjectToHandleCalOrCalStruct(calNow);
clear 'calNow';

%% Put them on common fitting basis, so that we are comparing the underlying
% data and not how it happened to be fit.
%
% Linear model basis
defaultNPrimaryBases = calNowOBJ.get('nPrimaryBases');
nPrimaryBases = input(sprintf('\nEnter number of primary bases [%d]: ',defaultNPrimaryBases));
if (isempty(nPrimaryBases))
    nPrimaryBases = defaultNPrimaryBases;
end
calThenOBJ.set('nPrimaryBases', nPrimaryBases);     % calThen.nPrimaryBases = nPrimaryBases;
calNowOBJ.set('nPrimaryBases', nPrimaryBases);      % calNow.nPrimaryBases = nPrimaryBases;
CalibrateFitLinMod(calThenOBJ);
CalibrateFitLinMod(calNowOBJ);

% Gamma type
defaultFitType = calNowOBJ.get('gamma.fitType');    % calNow.describe.gamma.fitType;
fitType = input(sprintf('Enter gamma fit type [%s]: ',defaultFitType),'s');
if (isempty(fitType))
    fitType = defaultFitType;
end
calThenOBJ.set('gamma.fitType', fitType);          % calThen.describe.gamma.fitType = fitType;
calNowOBJ.set('gamma.fitType', fitType);          % calNow.describe.gamma.fitType = fitType;
CalibrateFitGamma(calThenOBJ);
CalibrateFitGamma(calNowOBJ);

%% Say what we're doing
fprintf('\nComparing calibrations:\n');
fprintf('\t%s, %d, %s\n',thenFileName,thenIndex,calThenOBJ.get('date'));  % calThen.describe.date);
fprintf('\t%s, %d, %s\n',nowFileName,nowIndex,calNowOBJ.get('date'));    % calNow.describe.date);

%% Plot spectral power distributions.
%
% Plot as one plot if 3 or fewer primaries.
% Otherwise separate main measurements from what
% are probably the linear model correction terms.

then_S          = calThenOBJ.get('S');
then_P_device   = calThenOBJ.get('P_device');
then_P_ambient  = calThenOBJ.get('P_ambient');
then_gammaInput = calThenOBJ.get('gammaInput');
then_gammaTable = calThenOBJ.get('gammaTable');

now_S           = calNowOBJ.get('S');
now_P_device    = calNowOBJ.get('P_device');
now_P_ambient   = calNowOBJ.get('P_ambient');
now_gammaInput =  calNowOBJ.get('gammaInput');
now_gammaTable =  calNowOBJ.get('gammaTable');

now_nDevices    = calNowOBJ.get('nDevices');

if (size(now_gammaTable,2) <= now_nDevices) 
    figure; clf; hold on
    plot(SToWls(then_S), then_P_device, 'r');
    plot(SToWls(now_S),  now_P_device, 'g-');
    xlabel('Wavelength (nm)');
    ylabel('Power');
    title('Primaries');
else
    figure; clf;
    subplot(1,2,1); hold on
    plot(SToWls(then_S), then_P_device(:,1:now_nDevices),'r');
    plot(SToWls(now_S), now_P_device(:,1:now_nDevices),'g-');
    xlabel('Wavelength (nm)');
    ylabel('Power');
    title('Primaries');
    subplot(1,2,2); hold on
    plot(SToWls(then_S), then_P_device(:,now_nDevices+1:end),'r');
    plot(SToWls(now_S), now_P_device(:,now_nDevices+1:end),'g-');
    xlabel('Wavelength (nm)');
    ylabel('Power');
    title('Primaries (high order)');
end

%% Plot ambient
figure; clf; hold on
plot(SToWls(then_S), then_P_ambient,'r');
plot(SToWls(now_S), now_P_ambient,'g-');
xlabel('Wavelength (nm)');
ylabel('Power');
title('Ambient');

%% Explicitly compute and report ratio of R, G, and B full on spectra
rRatio = then_P_device(:,1)\now_P_device(:,1);
gRatio = then_P_device(:,2)\now_P_device(:,2);
bRatio = then_P_device(:,3)\now_P_device(:,3);
fprintf('Phosphor intensity ratios (now/then): %0.3g, %0.3g, %0.3g\n', ...
	rRatio,gRatio,bRatio);

%% Plot gamma functions
%
% Plot as one plot if 3 or fewer primaries.
% Otherwise separate main measurements from what
% are probably the linear model correction terms.
if (size(now_gammaTable,2) <= now_nDevices)
    figure; clf; hold on 
    plot(then_gammaInput, then_gammaTable,'r');
    plot(now_gammaInput, now_gammaTable,'g-');
    xlabel('Input');
    ylabel('Output');
    title('Gamma');
    ylim([0 1.2]);
else
    figure; clf;
    subplot(1,2,1); hold on
    plot(then_gammaInput, then_gammaTable(:,1:now_nDevices),'r');
    plot(now_gammaInput, now_gammaTable(:,1:now_nDevices),'g-');
    xlabel('Input');
    ylabel('Output');
    title('Gamma');
    ylim([0 1.2]);
    subplot(1,2,2); hold on
    plot(then_gammaInput,then_gammaTable(:,now_nDevices+1:end),'r');
    plot(now_gammaInput,now_gammaTable(:,now_nDevices+1:end),'g-');
    xlabel('Input');
    ylabel('Output');
    title('Gamma (high order)');
    ylim([-1.2 1.2]);
end

%% Let's print some luminance information
load T_xyzJuddVos;
T_xyz = SplineCmf(S_xyzJuddVos,683*T_xyzJuddVos,then_S);
S_xyz = then_S;
lumsThen = T_xyz(2,:)*then_P_device;
maxLumThen = sum(lumsThen(1:now_nDevices));
lumsNow = T_xyz(2,:)*now_P_device;
maxLumNow = sum(lumsNow(1:now_nDevices));
fprintf('Maximum luminance summing primaries: then %0.3g; now %0.3g\n',maxLumThen,maxLumNow);
minLumThen = T_xyz(2,:)*then_P_ambient;
minLumNow = T_xyz(2,:)*now_P_ambient;
fprintf('Minimum luminance: then %0.3g; now %0.3g\n',minLumThen,minLumNow);

%% Get max lum using calibration routines
SetSensorColorSpace(calThenOBJ,T_xyz,S_xyz);
SetSensorColorSpace(calNowOBJ,T_xyz,S_xyz);
maxXYZThen1 = SettingsToSensor(calThenOBJ,[1 1 1]');
maxXYZThen2 = SettingsToSensorAcc(calThenOBJ,[1 1 1]');
maxXYZNow1 = SettingsToSensor(calNowOBJ,[1 1 1]');
maxXYZNow2 = SettingsToSensorAcc(calNowOBJ,[1 1 1]');
fprintf('Maximum luminance SettingsToSensor: then %0.3g; now %0.3g\n',maxXYZThen1(2),maxXYZNow1(2));
fprintf('Maximum luminance SettingsToSensorAcc: then %0.3g; now %0.3g\n',maxXYZThen2(2),maxXYZNow2(2));

%% Plot new and old white point and channel chromaticities
figure; clf; hold on
maxxyYThen = XYZToxyY(maxXYZThen1);
maxxyYNow = XYZToxyY(maxXYZNow1);
plot(maxxyYThen(1),maxxyYThen(2),'ro','MarkerFaceColor','r','MarkerSize',10);
plot(maxxyYNow(1),maxxyYNow(2),'go','MarkerFaceColor','g','MarkerSize',10);

redxyYThen = XYZToxyY(SettingsToSensor(calThenOBJ,[1 0 0]'));
greenxyYThen = XYZToxyY(SettingsToSensor(calThenOBJ,[0 1 0]'));
bluexyYThen = XYZToxyY(SettingsToSensor(calThenOBJ,[0 0 1]'));
redxyYNow = XYZToxyY(SettingsToSensor(calNowOBJ,[1 0 0]'));
greenxyYNow = XYZToxyY(SettingsToSensor(calNowOBJ,[0 1 0]'));
bluexyYNow = XYZToxyY(SettingsToSensor(calNowOBJ,[0 0 1]'));
plot(redxyYThen(1),redxyYThen(2),'ro','MarkerFaceColor','r','MarkerSize',10);
plot(redxyYNow(1),redxyYNow(2),'go','MarkerFaceColor','g','MarkerSize',10);
plot(greenxyYThen(1),greenxyYThen(2),'ro','MarkerFaceColor','r','MarkerSize',10);
plot(greenxyYNow(1),greenxyYNow(2),'go','MarkerFaceColor','g','MarkerSize',10);
plot(bluexyYThen(1),bluexyYThen(2),'ro','MarkerFaceColor','r','MarkerSize',10);
plot(bluexyYNow(1),bluexyYNow(2),'go','MarkerFaceColor','g','MarkerSize',10);
axis('square');
axis([0.0 0.8 0.0 0.8]);
xlabel('x chromaticity');
ylabel('y chromaticity');
