% t_rgcEccData
% 
% Plot and fit data captured from plots in Coner & Kaplan (1995) and Dacey
% (2004) on RGC RF size as a function of eccentricity.
% 
% This data is used to scale ISETBIO RGC RF sizes as a function of
% eccentricity.
% 
% Data was digitized using 'digitize2' and ginput:
% https://www.mathworks.com/matlabcentral/fileexchange/928-digitize2-m
% 
% Croner & Kaplan, 1995, Figure 4:
% https://pdfs.semanticscholar.org/be5f/d167a456b4bce52e7edb7bd187616e79adc6.pdf
% 
% Dacey, 2004, Figure 2A:
% https://pdfs.semanticscholar.org/5f15/a3de07ccdbf2ef3763da262ede1c876f6b6a.pdf

%%
movingAverageFlag = true;

%% Dacey Parasol Data

load([isetbioRootPath '/isettools/data/rgc/parasolData.mat'])
vcNewGraphWin([],'wide'); 
subplot(1,3,1);
scatter(parasolData(:,1),parasolData(:,2))

% Linear regression
parasolFit = ([ones(size(parasolData,1),1) parasolData(:,1)]\parasolData(:,2));
% Fit regression
hold on; plot(.1:.1:15,(.1:.1:15).*parasolFit(2)+parasolFit(1));

% Moving average bin
if movingAverageFlag
xRange = 1;
for xRangeInd = 1:floor(max(parasolData(:,1))/xRange)
    xRangeVal(xRangeInd) = xRangeInd*xRange;
    xRangePts = find(abs(parasolData(:,1)-xRangeVal(xRangeInd))<xRange);
    if ~isempty(xRangePts)
    sizeAvg(xRangeInd) = mean(parasolData(xRangePts,2));
    end
    clear xRangePts
end
hold on;plot(xRangeVal,sizeAvg,'b')
end

title('Parasol DF Size over Eccentricity, Dacey 2004');
xlabel('Eccentricity (mm)'); ylabel(sprintf('Dendritic field size (\\mum)'));
grid on;
set(gca,'fontsize',14);

% Format plot
legend('Data','Fit','Binned Average');
axis([0 18 0 450]);


%% Croner and Kaplan parasol data

% This will not agree exactly with Dacey data which is dendritic field
% measurements

% figure; scatter([5 15 25]*.3,[.1 .18 .23]*.3)
% 
% midgetFit2 = ([ones(size([5 15 25]')) .3*[5 15 25]'])\(.3*[.1 .18 .23]')
% x = .1:.1:25; y = midgetFit2(1) + midgetFit2(2)*x;
% hold on; scatter(x,y)

% clear 
load([isetbioRootPath '/isettools/data/rgc/croner_kaplan_parasol_rgc.mat'])

% Microns to degrees
d1 = .2253*d1; 

% Radius (mm) to diameter (um)
d1(:,2) = 2*d1(:,2)*1000;

% figure; 
hold on;
scatter(d1(:,1),d1(:,2),'r')

% Linear regression
parasolFit = ([ones(size(d1,1),1) d1(:,1)]\d1(:,2));

% Plot regression
hold on; plot(.1:.1:8,(.1:.1:8).*parasolFit(2)+parasolFit(1));


% % Moving average bin
% if movingAverageFlag
% xRange = 5;
% for xRangeInd = 1:floor(max(d1(:,1))/xRange)
%     xRangeVal(xRangeInd) = xRangeInd*xRange;
%     xRangePts = find(abs(d1(:,1)-xRangeVal(xRangeInd))<xRange);
%     if ~isempty(xRangePts)
%     sizeAvg(xRangeInd) = mean(d1(xRangePts,2));
%     end
%     clear xRangePts
% end
% hold on;plot(xRangeVal,sizeAvg,'b')
% end

% Format plot
title(sprintf('Parasol DF Size over Eccentricity,\nDacey 2004 (blue) + Croner & Kaplan 1995 (red)'));
xlabel('Eccentricity (mm)'); ylabel(sprintf('Dendritic field size (\\mum)'));
grid on;
set(gca,'fontsize',14);

legend('Data','Fit','Binned Average');
% axis([0 18 0 450]);

%% Combine parasol data from Dacey and Croner & Kaplan

% parasolX = [parasolData(:,1); d1(:,1)];
% parasolY = [parasolData(:,2); d1(:,2)];
% figure; scatter(parasolX,parasolY);

%% Dacey Midget data

load([isetbioRootPath '/isettools/data/rgc/midgetData.mat'])
% figure; 
subplot(1,3,2);
scatter(midgetData(:,1),midgetData(:,2))

% Linear regression
midgetFit = ([ones(size(midgetData,1),1) midgetData(:,1)]\midgetData(:,2));
% Plot regression
hold on; plot(.1:.1:15,(.1:.1:15).*midgetFit(2)+midgetFit(1));

% Moving average bin
if movingAverageFlag
clear xRange xRangeVal sizeAvg
xRange = 1;
for xRangeInd = 1:floor(max(midgetData(:,1))/xRange)
    xRangeVal(xRangeInd) = xRangeInd*xRange;
    xRangePts = find(abs(midgetData(:,1)-xRangeVal(xRangeInd))<xRange);
    if ~isempty(xRangePts)
    sizeAvg(xRangeInd) = mean(midgetData(xRangePts,2));
    end
    clear xRangePts
end
hold on;plot(xRangeVal,sizeAvg,'b')
end

% Format plot
title('Midget DF Size over Eccentricity, Dacey 2004');
xlabel('Eccentricity (mm)'); ylabel(sprintf('Dendritic field size (\\mum)'));
grid on;
set(gca,'fontsize',14);

legend('Data','Fit','Binned Average');
axis([0 18 0 450]);


%% Dacey SBC data

load([isetbioRootPath '/isettools/data/rgc/sbcData.mat'])
% figure; 
subplot(1,3,3);
scatter(sbcData(:,1),-sbcData(:,2))

% Linear regression
sbcFit = ([-ones(size(sbcData,1),1) -sbcData(:,1)]\sbcData(:,2));
% Fit regression
hold on; plot(.1:.1:15,(.1:.1:15).*sbcFit(2)+sbcFit(1));
sbcx = sbcData(:,1); sbcy = -sbcData(:,2);

% Moving average bin
if movingAverageFlag
clear xRange xRangeVal sizeAvg
xRange = 1;
for xRangeInd = 1:floor(max(sbcData(:,1))/xRange)
    xRangeVal(xRangeInd) = xRangeInd*xRange;
    xRangePts = find(abs(sbcData(:,1)-xRangeVal(xRangeInd))<xRange);
    if ~isempty(xRangePts)
    sizeAvg(xRangeInd) = mean(-sbcData(xRangePts,2));
    end
    clear xRangePts
end
hold on;plot(xRangeVal,sizeAvg,'b')
end

% Format plot
title('SBC DF Size over Eccentricity, Dacey 2004');
xlabel('Eccentricity (mm)'); ylabel(sprintf('Dendritic field size (\\mum)'));
grid on;
set(gca,'fontsize',14);
legend('Data','Fit','Binned Average');
axis([0 18 0 450]);

%% Croner and Kaplan midget RF data
% Rough fit
plotFlag = false;
% This will not agree exactly with Dacey data which is dendritic field
% measurements

if plotFlag
deg2mm = 0.3;
figure; scatter([2.5 7.5 15 25 35]*deg2mm,[.03 .05 .07 .09 .15]*deg2mm)

% midgetFit2 = ([ones(size([2.5 7.5 15 25 35]')) [.03 .05 .07 .09 .15]'])\[2.5 7.5 15 25 35]';

% Linear regression
midgetFit2 = ([ones(size([2.5 7.5 15 25 35]')) deg2mm*[2.5 7.5 15 25 35]'])\(deg2mm*[.03 .05 .07 .09 .15]')
% Plot regression
x = .1:.1:35*.3; y = midgetFit2(1) + midgetFit2(2)*x;
hold on; scatter(x,y);
end


%% Watson RGC Formula - 2014


%% Get lookup table for how large RF is 
% There is an asymmetry in the size of RGC RFs over the retina
szCols = 128; fovRows = 90; fovCols = 90; scaleFactor = 1;
cellType = 'Midget';

% Top page 7, right column: "Typically we are concerned with the spacing
% within one class, in which case density is halved and the spacings should
% be multiplied by sqrt(2)."
% rgcDiameterLUT = scaleFactor*sqrt(2)*watsonRGCSpacing(szCols,szCols,fovRows)';
[rgcDiameterLUT, radDeg, rgc1d] = watsonRGCSpacing(szCols,szCols,fovRows);

figure; plot(radDeg,rgc1d); 
 xlabel('Eccentricity (degrees)'); ylabel('RF Size (degrees)'); grid on;
legend('Nasal','Superior','Temporal','Inferior');
 
title(sprintf('Human %s RGC RF Size (degrees)',cellType));
% figure; degStart = -fovCols/2; degEnd = fovCols/2;
% degarr = [degStart: (degEnd-degStart)/szCols : degEnd];
% contourf(degarr,degarr,rgcDiameterLUT,[0:max(rgcDiameterLUT(:))/20:max(rgcDiameterLUT(:))] ); axis square
% % surfc(degarr,degarr,rgcDiameterLUT); shading flat; %[0:max(rgcDiameterLUT(:))/20:max(rgcDiameterLUT(:))] ); % axis square
% title(sprintf('Human %s RGC RF Size (degrees)',cellType)); colorbar; 

%%
% figure; hold on;
% plot(rgcDiameterLUT(65:-1:1,65));
% plot(rgcDiameterLUT(65:end,65));
% plot(rgcDiameterLUT(65,65:-1:1));
% plot(rgcDiameterLUT(65,65:end));
% grid on;

%% Bipolar 

% Parasol and SBC

vcNewGraphWin;

% Parasol & SBC
subplot(1,2,1);
bipolarsPerRGC = (2 + (3/10)*(radDeg));

plot(radDeg,rgc1d(1,:).*bipolarsPerRGC/sqrt(2)); grid on;
axis([0 100 0 16]);
title('Diffuse Bipolar: Parasol and SBC');
 xlabel('Eccentricity (degrees)'); ylabel('RF Size (degrees)'); grid on;
 
% Midget
subplot(1,2,2);

bipolarsPerRGC = (1 + (2/10)*(radDeg));
plot(radDeg,rgc1d(1,:).*bipolarsPerRGC/sqrt(2)); grid on;
axis([0 100 0 16]);
title('Midget Bipolar');
 xlabel('Eccentricity (degrees)'); ylabel('RF Size (degrees)'); grid on;
% scatter(midgetData(:,1),midgetData(:,2))
