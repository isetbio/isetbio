%
% PAL_SDT_DPtoPCcomparison_Demo is a script for comparing the proportion 
% correct calculated for various d' (d-prime) for a 1AFC, standard 
% 2AFC, 2AFC same-different and 2AFC match-to-sample task, the last three 
% tasks under the Differencing model and all tasks assuming an unbiased 
% observer
%
% Syntax: PAL_SDT_DPtoPCcomparison_Demo
%
% asks for a vector of d' and outputs a table and a graph of proportion 
% correct against d' for a 1AFC, standard 2AFC, 2AFC same-different and 
% 2AFC match-to-sample task
%
% Example:
% 
% PAL_SDT_DPtoPCcomparison_Demo
% Enter a vector of d-prime values [0:.25:4]
%
% returns:
% 
%             -------------Proportion correct-----------
%     d-prime   1AFC      2AFC  2AFCsameDiff 2AFCmatchSamp
%          0    0.5000    0.5000    0.5000    0.5000
%     0.2500    0.5497    0.5702    0.5049    0.5057
%     0.5000    0.5987    0.6382    0.5195    0.5223
%     0.7500    0.6462    0.7021    0.5427    0.5486
%     1.0000    0.6915    0.7602    0.5733    0.5825
%     1.2500    0.7340    0.8116    0.6095    0.6216
%     1.5000    0.7734    0.8556    0.6495    0.6635
%     1.7500    0.8092    0.8920    0.6912    0.7058
%     2.0000    0.8413    0.9214    0.7330    0.7468
%     2.2500    0.8697    0.9442    0.7734    0.7850
%     2.5000    0.8944    0.9615    0.8110    0.8196
%     2.7500    0.9154    0.9741    0.8452    0.8501
%     3.0000    0.9332    0.9831    0.8753    0.8765
%     3.2500    0.9479    0.9892    0.9013    0.8989
%     3.5000    0.9599    0.9933    0.9231    0.9178
%     3.7500    0.9696    0.9960    0.9411    0.9336
%     4.0000    0.9772    0.9977    0.9555    0.9467
%
%FK (September 2009)

clear all;

if exist('OCTAVE_VERSION');
    fprintf('\nUnder Octave, Figure does not render exactly as intended. Visit\n');
    fprintf('www.palamedestoolbox.org/demosfiguregallery.html to see figure\n');
    fprintf('as intended.\n\n');
end

dP=input('Enter a vector of d-prime values ');


% Calculate PCs for a 1AFC task with optimum criterion=0
pHF=PAL_SDT_1AFC_DPtoPHF(dP,0);
pc1AFC=(pHF(:,1)+1.0-pHF(:,2))./2;
pc1AFC=pc1AFC';

% Calculate PCs for a 2AFC task
pc2AFC=PAL_SDT_MAFC_DPtoPC(dP,2);

% Calculate PCs for a 2AFC Same Different task
pc2AFCsameDiff=PAL_SDT_2AFCsameDiff_DPtoPC(dP);

% Calculate PCs for a 2AFC Match-to_Sample task using a Differencing strategy 
pc2AFCmatchSample=PAL_SDT_2AFCmatchSample_DiffMod_DPtoPC(dP);

% Type out results
output=[dP;pc1AFC;pc2AFC;pc2AFCsameDiff;pc2AFCmatchSample];
output=output';
fprintf('\n');
disp('              -------------Proportion correct-----------');
disp('    d-prime   1AFC      2AFC  2AFCsameDiff 2AFCmatchSamp');
disp(output);

% Plot a graph of the results
figure('name','Pc versus d-prime for various tasks');
set(gca, 'fontsize',18);
plot(dP,pc2AFC, 'g-s','MarkerEdgeColor','k','MarkerFaceColor','k','Markersize',10,'LineWidth',2);
hold on;
plot(dP,pc1AFC, 'g-d','MarkerEdgeColor','k','MarkerFaceColor','k','Markersize',10,'LineWidth',2);
hold on;
plot(dP,pc2AFCmatchSample, 'g-x','MarkerEdgeColor','k','MarkerFaceColor','k','Markersize',10,'LineWidth',2);
hold on;
plot(dP,pc2AFCsameDiff, 'g-o','MarkerEdgeColor','k','MarkerFaceColor','k','Markersize',10,'LineWidth',2);


% Add some labeling
xlabel('d-prime','FontSize',18);
ylabel('Proportion correct','FontSize',18);
h = legend('2AFC','1AFC','2AFC match-to-sample','2AFC same-different',2);
set(h,'Interpreter','none','fontsize',16,'Location','SouthEast');