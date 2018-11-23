% Watson RGC Formula - 2014
%
% Description:
%    The retinal ganglion cell density estimates from 
%       Watson 2014: A formula for human retinal ganglion cell receptive
%       field density as a function of visual field location  
%       http://jov.arvojournals.org/article.aspx?articleid=2279458
%
% See Also:
%    t_rgcEccData
%

% History:
%    XX/XX/17  JRG/BW  ISETBIO Team, 2017
%    11/22/18  JNM     Formatting

%% Get lookup table for how large RF is 
% There is an asymmetry in the size of RGC RFs over the retina
szCols = 128;
fovRows = 90;
fovCols = 90;
scaleFactor = 1;
cellType = 'Midget';

% Top page 7, right column: "Typically we are concerned with the spacing
% within one class, in which case density is halved and the spacings should
% be multiplied by sqrt(2)."
% rgcDiameterLUT = scaleFactor * sqrt(2) ...
%    * watsonRGCSpacing(szCols, szCols, fovRows)';
% [rgcDiameterLUT, radDeg, rgc1d] = ...
%    watsonRGCSpacing(szCols, szCols, fovRows);
[rgcDiameterLUT, radDeg, rgc1d] = watsonRGCSpacing(fovRows);

%% Plots
arcMinPerDegree = 60;
convertDensityFactor = sqrt(2);
vcNewGraphWin;
cind = 'rbgk';
hold on;
for k = 1:4
    plot(radDeg, convertDensityFactor * arcMinPerDegree * rgc1d(k, :), ...
        cind(k), 'linewidth', 2);
end
xlabel('Eccentricity (degrees)');
ylabel('RF Size (degrees)');
grid on;
legend('Temporal', 'Superior', 'Nasal', 'Inferior', 'location', 'nw')
title(sprintf('Human %s RGC RF Size (degrees)', cellType));
axis([0 10 0 6]);
set(gca, 'fontsize', 14);

% figure;
% degStart = -fovCols / 2;
% degEnd = fovCols / 2;
% degarr = [degStart: (degEnd - degStart) / szCols : degEnd];
% contourf(degarr, degarr, convertDensityFactor * rgcDiameterLUT, ...
%    [0:max(rgcDiameterLUT(:)) / 20:max(rgcDiameterLUT(:))] );
% axis square
% % surfc(degarr, degarr, rgcDiameterLUT);
% shading flat;
% %[0:max(rgcDiameterLUT(:)) / 20:max(rgcDiameterLUT(:))] );
% % axis square
% title(sprintf('Human %s RGC RF Size (degrees)', cellType));
% colorbar;

%%