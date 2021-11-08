%% Enchroma transmittance
%
% The Enchroma glasses are an interesting attempt to help anomalous
% trichromats experience more color.  Class projects with these glasses can
% be informative.
%
% Data were taken from Figure 6 in this paper
%
% https://www.semanticscholar.org/paper/Assessment-of-Enchroma-Filter-for-Correcting-Color-Almutairi-Kundart/b8cbf93781f85db40bc404755111ef324e69836d
%
% Assessment of Enchroma Filter for Correcting Color Vision Deficiency
% Almutairi et al.
% 

%% Data grabbed with grabit.  Sorry for the bad names.

load('EnchromaInput','EnchromaInput');
load('EnchromaThroughLens','EnchromGrabThroughLens');

%% Interpolate - default method is linear.  Seems OK.

wave = (350:10:1000);
enchromIn = interp1(EnchromaInput(:,1), EnchromaInput(:,2),wave);
ieNewGraphWin; 
plot(EnchromaInput(:,1),EnchromaInput(:,2));

enchromOut = interp1(EnchromGrabThroughLens(:,1), EnchromGrabThroughLens(:,2),wave);
ieNewGraphWin; 
plot(EnchromGrabThroughLens(:,1),EnchromGrabThroughLens(:,2));

%% Plot the transmittance

% Force the small bad values to within 0,1 range.
transmittance = enchromOut./enchromIn;
transmittance = min(transmittance,1);
transmittance = max(transmittance,0);

ieNewGraphWin;
plot(wave,transmittance,'-o');
grid on;
xlabel('Wave (nm)');
ylabel('Transmittance');

%% Save and check that it worked

fname = fullfile(isetRootPath,'data','lens','enchroma.mat');
ieSaveSpectralFile(wave(:),transmittance(:),'Estimated Enchroma transmittance',fname);

tmp = ieReadSpectra(fname,wave);
ieNewGraphWin;
plot(wave,tmp);

%%
