%%coneDensitySong2011MakeData
%
% Description:
%    Script to make data files for cone density from paper:
%      Song, H., Chui, T. Y. P., Zhong, Z., Elsner, A. E., & Burns, S. A.
%      (2011). Variation of Cone Photoreceptor Packing Density with Retinal
%      Eccentricity and Age. Investigative Ophthalmology & Visual Science,
%      52(10), 7376-7384. http://doi.org/10.1167/iovs.11-7199
%
%    These data are from Table 1 of the paper, copied and pasted.
%    Order of the rows is Superior, Inferior, Nasal, Temporal  (SINT).
%    There are separate data for old and young subjects.
%
%    These numbers are substantially lower in the central fovea than the
%    Curcio numbers that are published (around 250,000 and here around
%    70,000).

% BW Vistasoft Team, 2016
%
% 08/16/17  dhb  Rename, more comments.

%% Initialize
ieInit;

%% Where the data will be saved
chdir(fullfile(isetbioRootPath,'data','datafiles','human','cones'));

%% Young subject data
%
% The commented out version has error bars.
%
% YoungRaw = [63.5 ± 2.3	52.7 ± 1.9	42.4 ± 2.0	34.0 ± 2.1	29.6 ± 1.9	22.1 ± 1.4 18.5 ± 1.2	15.7 ± 0.8	13.1 ± 0.8	11.7 ± 0.9	10.0 ± 0.8	8.7 ± 0.5;
%     62.8 ± 4.5	55.3 ± 1.4	44.6 ± 3.3	36.3 ± 2.8	31.4 ± 2.7	23.9 ± 2.1 19.4 ± 1.6	16.6 ± 1.7	12.8 ± 1.1	11.5 ± 0.9	10.2 ± 0.8	8.1 ± 0.7;
%     68.2 ± 5.4	59.7 ± 2.8	50.0 ± 2.9	43.7 ± 2.7	37.8 ± 2.2	29.1 ± 1.8 24.2 ± 1.2	19.1 ± 1.5	16.8 ± 1.5	14.5 ± 1.4	11.9 ± 0.9	10.4 ± 0.6;
%     75.2 ± 7.5	59.2 ± 4.0	50.5 ± 2.4	41.2 ± 2.0	37.3 ± 1.7	28.1 ± 1.3 24.1 ± 1.5	19.9 ± 1.4	16.3 ± 0.9	13.2 ± 1.0	11.5 ± 0.5	9.7 ± 0.7];
eccentricityMm = [0.18, 0.27, 0.36, 0.45, 0.54, 0.72, 0.90, 1.08, 1.35, 1.62, 1.89, 2.16]; 
YoungRaw = ...
 [ 63.5000   52.7000   42.4000   34.0000   29.6000   22.1000   18.5000   15.7000   13.1000   11.7000   10.0000    8.7000;
   62.8000   55.3000   44.6000   36.3000   31.4000   23.9000   19.4000   16.6000   12.8000   11.5000   10.2000    8.1000;
   68.2000   59.7000   50.0000   43.7000   37.8000   29.1000   24.2000   19.1000   16.8000   14.5000   11.9000   10.4000;
   75.2000   59.2000   50.5000   41.2000   37.3000   28.1000   24.1000   19.9000   16.3000   13.2000   11.5000    9.7000];
vcNewGraphWin; plot(eccentricityMm,YoungRaw)
grid on; xlabel('mm'); ylabel('Density 10^3 cones/mm^2');
title('Song et al. estimates, young subjects');

%% Old subject data
%
% The commented out version has error bars.
%
% OldRaw = [50.2 ± 2.1	46.5 ± 2.4	39.9 ± 3.1	33.7 ± 1.9	29.4 ± 1.7	23.2 ± 1.6 19.0 ± 0.9	16.4 ± 1.4	12.8 ± 0.8	11.1 ± 0.7	10.0 ± 0.6	9.1 ± 0.8;
%     50.1 ± 1.7	43.6 ± 2.2	39.5 ± 1.5	31.2 ± 1.6	28.8 ± 1.7	22.8 ± 1.6 18.4 ± 1.2	15.2 ± 0.8	11.3 ± 0.7	11.1 ± 0.8	8.8 ± 0.6	8.3 ± 0.7;
%     52.6 ± 4.4	48.3 ± 2.7	43.0 ± 1.8	35.6 ± 2.3	33.1 ± 2.4	27.5 ± 1. 22.5 ± 1.1	20.9 ± 1.2	17.6 ± 1.3	15.4 ± 1.4	12.4 ± 0.8	13.0 ± 1.2;
%     46.6 ± 1.9	40.7 ± 1.6	39.0 ± 3.0	35.6 ± 2.8	33.6 ± 2.4	25.8 ± 1.9 22.0 ± 1.5	18.3 ± 1.2	14.9 ± 0.8	13.3 ± 0.6	11.0 ± 0.7	9.0 ± 0.9 ]
OldRaw = ...
  [50.2000   46.5000   39.9000   33.7000   29.4000   23.2000   19.0000   16.4000   12.8000   11.1000   10.0000    9.1000;
   50.1000   43.6000   39.5000   31.2000   28.8000   22.8000   18.4000   15.2000   11.3000   11.1000    8.8000    8.3000;
   52.6000   48.3000   43.0000   35.6000   33.1000   27.5000   22.5000   20.9000   17.6000   15.4000   12.4000   13.0000;
   46.6000   40.7000   39.0000   35.6000   33.6000   25.8000   22.0000   18.3000   14.9000   13.3000   11.0000    9.0000];

%% Plot old and young
vcNewGraphWin; plot(eccentricityMm,OldRaw)
grid on; xlabel('mm'); ylabel('Density 10^3 cones/mm^2');
title('Song et al. estimates, old subjects');

%% Convert to cones/mm2
Old   = OldRaw*10^3;
Young = YoungRaw*10^3;

%% Build and save cone density structure for young eyes
superior.eccMM = eccentricityMm;
superior.units = 'cones/mm2';
superior.density = Young(1,:);

inferior.eccMM = eccentricityMm;
inferior.units = 'cones/mm2';
inferior.density = Young(2,:);

nasal.eccMM = eccentricityMm;
nasal.units = 'cones/mm2';
nasal.density = Young(3,:);

temporal.eccMM = eccentricityMm;
temporal.units = 'cones/mm2';
temporal.density = Young(4,:);

description = 'From Table 1 in Variation of Cone Photoreceptor Packing Density with Retinal Eccentricity and Age, Song et al. IOVS; young eyes.';

save('coneDensitySong2011Young','superior','inferior','nasal','temporal','description');
dYoung = load('coneDensitySong2011Young');

%% Start a comparison plot
vcNewGraphWin; 
plot(inferior.eccMM(3:end),inferior.density(3:end),'-ro');
grid on;

%% Build and save cone density structure for old eyes
superior.eccMM = eccentricityMm;
superior.units = 'cones/mm2';
superior.density = Old(1,:);

inferior.eccMM = eccentricityMm;
inferior.units = 'cones/mm2';
inferior.density = Old(2,:);

nasal.eccMM = eccentricityMm;
nasal.units = 'cones/mm2';
nasal.density = Old(3,:);

temporal.eccMM = eccentricityMm;
temporal.units = 'cones/mm2';
temporal.density = Old(4,:);

description = 'From Table 1 in Variation of Cone Photoreceptor Packing Density with Retinal Eccentricity and Age, Song et al. IOVS; old eyes.';

save('coneDensitySong2011Old','superior','inferior','nasal','temporal','description');
dOld = load('coneDensitySong2011Old');

%% Add to comparison plot
hold on
plot(dOld.inferior.eccMM,dOld.inferior.density,'-go');
grid on;

%% Add in Curcio data to plot
dCurcio = load('coneDensityCurcio1990');
hold on
plot(dCurcio.inferior.eccMM,dCurcio.inferior.density,'-bo');
grid on;
legend({'young','old','curcio'})

%% Get Curcio data via coneDensityReadData and make sure it matches what is in the file
%
% Also shows how to pass eccentricity in mm to coneDensityReadData.
eccentricityMm = linspace(0,20,100);
curcio = zeros(size(eccentricityMm));
for ii=1:length(eccentricityMm)
    curcio(ii) = coneDensityReadData('eccentricity',eccentricityMm(ii),'angle',270,'whichEye','left','eccentricityUnits','mm');
end
vcNewGraphWin; hold on
plot(dCurcio.inferior.eccMM,dCurcio.inferior.density,'-bo');
plot(eccentricityMm,curcio,'gx');
title('Curcio data two ways');


%% Get Song data via coneDensityReadData and make sure it matches what is in the file
%
% Also shows how to pass eccentricity in microns to coneDensityReadData, and angle in radians
eccentricityMm = linspace(0,4,20);
songOld = zeros(size(eccentricityMm));
songYoung = zeros(size(eccentricityMm));
for ii=1:length(eccentricityMm)
    songOld(ii) = coneDensityReadData('eccentricity',1e3*eccentricityMm(ii),'angle',270,'whichEye','left','coneDensitySource','Song2011Old','eccentricityUnits','um');
    songYoung(ii) = coneDensityReadData('eccentricity',1e-3*eccentricityMm(ii),'angle',3*pi/2,'whichEye','left','coneDensitySource','Song2011Young','angleUnits','rad');
end
vcNewGraphWin; hold on
plot(dOld.inferior.eccMM,dOld.inferior.density,'-ro');
plot(eccentricityMm,songOld,'gx');
plot(dYoung.inferior.eccMM,dYoung.inferior.density,'-ko');
plot(eccentricityMm,songYoung,'bx');
title('Song data two ways');

