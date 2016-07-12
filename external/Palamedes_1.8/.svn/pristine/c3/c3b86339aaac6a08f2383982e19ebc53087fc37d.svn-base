%
% PAL_SDT_PCtoDPcomparison_Demo is a script for comparing d's 
% calculated for various proportion corrects for a 1AFC, standard 
% 2AFC, 2AFC same-different and 2AFC match-to-sample task, the last three 
% tasks under the Differencing model and all tasks assuming an unbiased 
% observer
%
% Syntax: PAL_SDT_PCtoDPcomparison_Demo
%
% asks for a vector of proportion correct and outputs a table and graph
% of d' against proportion correct for a for a 1AFC, standard 2AFC, 2AFC 
% same-different and 2AFC match-to-sample task
%
% Example:
%
% PAL_SDT_PCtoDPcomparison_Demo
%
% Enter a vector of proportion correct values [.5:.025:.975]
% 
%               -------------------- d-prime --------------
%     p Corr    1AFC      2AFC   2AFCsameDiff 2AFCmatchSamp
%     0.5000         0   -0.0000         0   -0.0000
%     0.5250    0.1254    0.0887    0.5680    0.5297
%     0.5500    0.2513    0.1777    0.8146    0.7614
%     0.5750    0.3782    0.2675    1.0124    0.9486
%     0.6000    0.5067    0.3583    1.1872    1.1152
%     0.6250    0.6373    0.4506    1.3490    1.2708
%     0.6500    0.7706    0.5449    1.5032    1.4206
%     0.6750    0.9075    0.6417    1.6535    1.5679
%     0.7000    1.0488    0.7416    1.8022    1.7153
%     0.7250    1.1955    0.8453    1.9515    1.8652
%     0.7500    1.3490    0.9539    2.1036    2.0200
%     0.7750    1.5108    1.0684    2.2605    2.1822
%     0.8000    1.6832    1.1902    2.4246    2.3549
%     0.8250    1.8692    1.3217    2.5990    2.5420
%     0.8500    2.0729    1.4657    2.7879    2.7493
%     0.8750    2.3007    1.6269    2.9972    2.9852
%     0.9000    2.5631    1.8124    3.2368    3.2631
%     0.9250    2.8791    2.0358    3.5243    3.6079
%     0.9500    3.2897    2.3262    3.8976    4.0724
%     0.9750    3.9199    2.7718    4.4730    4.8142
%
%FK (September 2009)

clear all;

if exist('OCTAVE_VERSION');
    fprintf('\nUnder Octave, Figure does not render exactly as intended. Visit\n');
    fprintf('www.palamedestoolbox.org/demosfiguregallery.html to see figure\n');
    fprintf('as intended.\n\n');
end

pC=input('Enter a vector of proportion correct values ');


% Calculate pHit and pFA for 1AFC task with optimum criterion (C=0)
pHF(:,1)=pC;
pHF(:,2)=1-pC;

% Calculate DPs for basic 1AFC task
dp1AFC=PAL_SDT_1AFC_PHFtoDP(pHF);
dp1AFC=dp1AFC';

% Calculate DPs for a 2AFC task
dp2AFC=PAL_SDT_MAFC_PCtoDP(pC,2);

% Calculate DPs for a 2AFC Same Different task
dp2AFCsameDiff=PAL_SDT_2AFCsameDiff_PCtoDP(pC);

% Calculate DPs for a 2AFC Match-to-Sample task assuming a Differencing
% strategy
dp2AFCmatchSample=PAL_SDT_2AFCmatchSample_DiffMod_PCtoDP(pC);

% Type out results
output=[pC;dp1AFC;dp2AFC;dp2AFCsameDiff;dp2AFCmatchSample];
output=output';
fprintf('\n');
disp('              -------------------- d-prime --------------');
disp('    p Corr    1AFC      2AFC   2AFCsameDiff 2AFCmatchSamp');
disp(output);

figure('name','d-prime versus pC for various tasks');
set(gca, 'fontsize',18);
plot(pC,dp2AFC, 'g-s','MarkerEdgeColor','k','MarkerFaceColor','k','Markersize',10,'LineWidth',2);
hold on;
plot(pC,dp1AFC, 'g-d','MarkerEdgeColor','k','MarkerFaceColor','k','Markersize',10,'LineWidth',2);
hold on;
plot(pC,dp2AFCsameDiff, 'g-x','MarkerEdgeColor','k','MarkerFaceColor','k','Markersize',10,'LineWidth',2);
hold on;
plot(pC,dp2AFCmatchSample, 'g-o','MarkerEdgeColor','k','MarkerFaceColor','k','Markersize',10,'LineWidth',2);

% Add some labeling
xlabel('Proportion correct','FontSize',18);
ylabel('d-prime','FontSize',18);
h = legend('2AFC','1AFC','2AFC same-different','2AFC match-to-sample',2);
set(h,'Interpreter','none','fontsize',16,'Location','NorthWest');