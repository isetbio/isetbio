%
% PAL_SDT_DPtoPCacrossM_Demo is a script for comparing the effect of 
% the number of stimulus alternatives M on proportion correct for a 
% Standard, Match-to-Sample and Oddity task, the last two tasks assuming 
% a Differencing model, and all three tasks assuming an unbiased observer
%
% Syntax: PAL_SDT_DPtoPCacrossM_Demo
%
% asks for a single (scalar) d' and outputs a table and a graph of 
% proportion correct against the  number of alternatives M for a Standard, 
% Match-to-Sample and Oddity task.  Note that the output for the last two
% tasks will not be identical each time since the computations involve
% Monte Carlo simulation
%
% Example:
% 
% PAL_SDT_DPtoPCacrossM_Demo
%
% Enter a single d-prime value [2]
%
% returns:
% 
% 
%              ----Proportion correct-----
%        M     Standard  MatchSamp  Oddity
%     3.0000    0.8658    0.6449    0.6040
%     4.0000    0.8228    0.5807    0.5985
%     5.0000    0.7878    0.5335    0.5863
%     6.0000    0.7585    0.4958    0.5722
%     7.0000    0.7332    0.4651    0.5572
%     8.0000    0.7110    0.4391    0.5422
%
% FK (September 2009)

clear all;

if exist('OCTAVE_VERSION');
    fprintf('\nUnder Octave, Figure does not render exactly as intended. Visit\n');
    fprintf('www.palamedestoolbox.org/demosfiguregallery.html to see figure\n');
    fprintf('as intended.\n\n');
end

dP=input('Enter a single d-prime value ');
disp('Working. Please wait.');

% Create some vectors
M=[3 4 5 6 7 8];
pcMAFCstandard=zeros(1,6);
pcMAFCmatchSample=zeros(1,6);
pcMAFCoddity=zeros(1,6);

% Calculate PCs for MAFC Standard, Match-to-Sample and Oddity tasks 
for n=1:length(M)
pcMAFCstandard(n)=PAL_SDT_MAFC_DPtoPC(dP,M(n));
pcMAFCmatchSample(n)=PAL_SDT_MAFCmatchSample_DiffMod_DPtoPC(dP,M(n));
pcMAFCoddity(n)=PAL_SDT_MAFCoddity_DPtoPC(dP,M(n));
end

% Type out results
output=[M;pcMAFCstandard;pcMAFCmatchSample;pcMAFCoddity];
output=output';
fprintf('\n');
disp('             ----Proportion correct-----');
disp('       M     Standard  MatchSamp  Oddity');
disp(output);

% Plot a graph of the results
figure('name','Pc from d-prime as a function of M');
set(gca, 'fontsize',16);
plot(M,pcMAFCstandard, 'g-s','MarkerEdgeColor','k','MarkerFaceColor','k','Markersize',10,'LineWidth',2);
hold on;
plot(M,pcMAFCmatchSample, 'g-d','MarkerEdgeColor','k','MarkerFaceColor','k','Markersize',10,'LineWidth',2);
hold on;
plot(M,pcMAFCoddity, 'g-x','MarkerEdgeColor','k','MarkerFaceColor','k','Markersize',10,'LineWidth',2);
hold on;
set(gca,'xtick',M);

% Add some labeling
xlabel('Number of alternatives M','FontSize',18);
ylabel('Proportion correct','FontSize',18);
h = legend('Standard','Match-to-Sample','Oddity',2);
set(h,'Interpreter','none','fontsize',16,'Location','NorthEast');