%% Variation of Cone Photoreceptor Packing Density with Retinal Eccentricity and Age
% Hongxin Song; Toco Yuen Ping Chui; Zhangyi Zhong; Ann E. Elsner; Stephen A. Burns
% Investigative Ophthalmology & Visual Science September 2011, Vol.52, 7376-7384. doi:10.1167/iovs.11-7199
%
% These data are from Table 1, copied and pasted.
%
% Order of the rows is
%  Superior, Inferior, Nasal, Temporal  (SINT)
%
% These numbers are substantially lower in the central fovea than the
% Curcio numbers that are published (around 250,000 and here around
% 70,000).
%
% BW Vistasoft Team, 2016

%% Where the data will be saved
chdir(fullfile(isetbioRootPath,'/isettools/data/human/coneMosaic/'));

%% With error bars
%
% Young = [63.5 ± 2.3	52.7 ± 1.9	42.4 ± 2.0	34.0 ± 2.1	29.6 ± 1.9	22.1 ± 1.4 18.5 ± 1.2	15.7 ± 0.8	13.1 ± 0.8	11.7 ± 0.9	10.0 ± 0.8	8.7 ± 0.5;
%     62.8 ± 4.5	55.3 ± 1.4	44.6 ± 3.3	36.3 ± 2.8	31.4 ± 2.7	23.9 ± 2.1 19.4 ± 1.6	16.6 ± 1.7	12.8 ± 1.1	11.5 ± 0.9	10.2 ± 0.8	8.1 ± 0.7;
%     68.2 ± 5.4	59.7 ± 2.8	50.0 ± 2.9	43.7 ± 2.7	37.8 ± 2.2	29.1 ± 1.8 24.2 ± 1.2	19.1 ± 1.5	16.8 ± 1.5	14.5 ± 1.4	11.9 ± 0.9	10.4 ± 0.6;
%     75.2 ± 7.5	59.2 ± 4.0	50.5 ± 2.4	41.2 ± 2.0	37.3 ± 1.7	28.1 ± 1.3 24.1 ± 1.5	19.9 ± 1.4	16.3 ± 0.9	13.2 ± 1.0	11.5 ± 0.5	9.7 ± 0.7];

ecc = [0.18, 0.27, 0.36, 0.45, 0.54, 0.72, 0.90, 1.08, 1.35, 1.62, 1.89, 2.16]; % mm
Young = ...
 [ 63.5000   52.7000   42.4000   34.0000   29.6000   22.1000   18.5000   15.7000   13.1000   11.7000   10.0000    8.7000;
   62.8000   55.3000   44.6000   36.3000   31.4000   23.9000   19.4000   16.6000   12.8000   11.5000   10.2000    8.1000;
   68.2000   59.7000   50.0000   43.7000   37.8000   29.1000   24.2000   19.1000   16.8000   14.5000   11.9000   10.4000;
   75.2000   59.2000   50.5000   41.2000   37.3000   28.1000   24.1000   19.9000   16.3000   13.2000   11.5000    9.7000];
vcNewGraphWin; plot(ecc,Young)
grid on; xlabel('mm'); ylabel('Density 10^3 cones/mm^2');

%% With error bars
%
% Old = ...
%     [50.2 ± 2.1	46.5 ± 2.4	39.9 ± 3.1	33.7 ± 1.9	29.4 ± 1.7	23.2 ± 1.6 19.0 ± 0.9	16.4 ± 1.4	12.8 ± 0.8	11.1 ± 0.7	10.0 ± 0.6	9.1 ± 0.8;
%     50.1 ± 1.7	43.6 ± 2.2	39.5 ± 1.5	31.2 ± 1.6	28.8 ± 1.7	22.8 ± 1.6 18.4 ± 1.2	15.2 ± 0.8	11.3 ± 0.7	11.1 ± 0.8	8.8 ± 0.6	8.3 ± 0.7;
%     52.6 ± 4.4	48.3 ± 2.7	43.0 ± 1.8	35.6 ± 2.3	33.1 ± 2.4	27.5 ± 1. 22.5 ± 1.1	20.9 ± 1.2	17.6 ± 1.3	15.4 ± 1.4	12.4 ± 0.8	13.0 ± 1.2;
%     46.6 ± 1.9	40.7 ± 1.6	39.0 ± 3.0	35.6 ± 2.8	33.6 ± 2.4	25.8 ± 1.9 22.0 ± 1.5	18.3 ± 1.2	14.9 ± 0.8	13.3 ± 0.6	11.0 ± 0.7	9.0 ± 0.9 ]

Old = ...
  [50.2000   46.5000   39.9000   33.7000   29.4000   23.2000   19.0000   16.4000   12.8000   11.1000   10.0000    9.1000;
   50.1000   43.6000   39.5000   31.2000   28.8000   22.8000   18.4000   15.2000   11.3000   11.1000    8.8000    8.3000;
   52.6000   48.3000   43.0000   35.6000   33.1000   27.5000   22.5000   20.9000   17.6000   15.4000   12.4000   13.0000;
   46.6000   40.7000   39.0000   35.6000   33.6000   25.8000   22.0000   18.3000   14.9000   13.3000   11.0000    9.0000];


vcNewGraphWin; plot(ecc,Old)
grid on; xlabel('mm'); ylabel('Density 10^3 cones/mm^2');

%% Plot the Curcio data

% Convert to cones/mm2
Old   = Old*10^3;
Young = Young*10^3;

%% Build cone density structure for young eyes
d.superior.eccMM = ecc;
d.superior.units = 'cones/mm2';
d.superior.density = Young(1,:);

d.inferior.eccMM = ecc;
d.inferior.units = 'cones/mm2';
d.inferior.density = Young(2,:);

d.nasal.eccMM = ecc;
d.nasal.units = 'cones/mm2';
d.nasal.density = Young(3,:);

d.temporal.eccMM = ecc;
d.temporal.units = 'cones/mm2';
d.temporal.density = Young(4,:);

d.description = 'From Table 1 in Variation of Cone Photoreceptor Packing Density with Retinal Eccentricity and Age, Song et al. IOVS';

save('coneDensityYoung','d');

vcNewGraphWin; 
plot(d.inferior.eccMM(3:end),d.inferior.density(3:end),'-ro');
grid on;

%% Build cone density structure for old eyes
d.superior.eccMM = ecc;
d.superior.units = 'cones/mm2';
d.superior.density = Old(1,:);

d.inferior.eccMM = ecc;
d.inferior.units = 'cones/mm2';
d.inferior.density = Old(2,:);

d.nasal.eccMM = ecc;
d.nasal.units = 'cones/mm2';
d.nasal.density = Old(3,:);

d.temporal.eccMM = ecc;
d.temporal.units = 'cones/mm2';
d.temporal.density = Old(4,:);

d.description = 'From Table 1 in Variation of Cone Photoreceptor Packing Density with Retinal Eccentricity and Age, Song et al. IOVS';

save('coneDensityOld','d');

hold on
plot(d.inferior.eccMM(3:end),d.inferior.density(3:end),'-go');
grid on;

%%

d = load('coneDensity');
hold on
plot(d.inferior.eccMM(3:end),d.inferior.density(3:end),'-bo');
grid on;

legend({'young','old','curcio'})


pp = photoPigment;
ecc = 0:0.1:2;
curcio = zeros(size(ecc));
burns  = zeros(size(ecc));
for ii=1:length(ecc)
    curcio(ii) = pp.eccDensity('eccMM',ecc(ii));
    burns(ii)  = pp.eccDensity('eccMM',ecc(ii),'dataSet','burns young');
end

vcNewGraphWin;
plot(ecc,burns,'r-o',ecc,curcio,'g-x');
grid on;

