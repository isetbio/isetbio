%
% PAL_SDT_ROCML_Demo is a script for demonstrating the fitting of three 
% ROC curves derived from 1AFC rating-scale data, using a Maximum 
% Likelihood criterion. The raw data for the highest d-prime value is 
% taken from Table 5.1 in McNicol (2005) and for the middle d-prime value 
% from Table 3.5 in Macmillan & Creelman (2005). The corresponding output 
% graphs are similar to those in Fig. 5.1 of McNicol (2005) and Fig. 3.2 of 
% Macmillan & Creelman (2005).
% 
% Refs: 
% 
% Macmillan, N.A. & Creelman, C.D. (2005) Detection Theory: A User's 
% Guide: Lawrence Erlbaum Associates.
%
% McNicol, D. (2005) A Primer of Signal Detection Theory: Lawrence Erlbaum
% Associates.
%
% Syntax: PAL_SDT_ROCML_Demo
%
%
% 
% Introduced: Palamedes version 1.6.0 (FK & NP)
%

clear all;

message = 'Enter number of simulations to determine standard errors (e.g. 400) ';
Bse = input(message);
message = 'Enter number of simulations to determine goodness-of-fit (e.g. 400) ';
Bgof = input(message);
message = 'Enter number of simulations to determine ratio-of-SDs comparison (e.g. 400) ';
Brsd = input(message);

SDTF = @PAL_SDT_1AFC_PHFtoDP;
invSDTF = @PAL_SDT_1AFC_DPtoPHF;

options = PAL_minimize('options');
options.TolX = 1e-4;        %sacrifice precision for speed (default: 1e-6)
options.TolFun = 1e-4;      %sacrifice precision for speed (default: 1e-6)

% Create three matrices of raw numbers of 'hits' and 'false alarms', each 
% for five different ratings,in order to derive three ROC curves.  Each row
% contains a 'hit' and its corresponding 'false alarm'
NumHF_1 = [159 2; 41 3; 19 21; 37 80; 32 182];
NumHF_2 = [37 4; 25 18; 18 28; 11 21; 9 29];
NumHF_3 = [103 30; 92 48; 73 60; 54 117; 53 120];

% The following routines take the raw numbers and converts them into
% cumulative scores, and creates matrices of the total number of trials
% and the proportion of hits and false alarms
[cumNumHF_1 OutOfNum_1 pHF_1] = PAL_SDT_cumulateHF(NumHF_1);
[cumNumHF_2 OutOfNum_2 pHF_2] = PAL_SDT_cumulateHF(NumHF_2);
[cumNumHF_3 OutOfNum_3 pHF_3] = PAL_SDT_cumulateHF(NumHF_3);

% Estimate the best fitting parameters of d-prime and SD ratio for each ROC
% with initial guess for the SD ratio set to 0.75 or 1 (this is an optional
% argument that if not included defaults to 1)  
[dP_1 R_1 C_1 negLL exitflag] = PAL_SDT_ROCML_Fit(cumNumHF_1, OutOfNum_1,SDTF,invSDTF,'ratioSDvalue',0.75,'searchOptions',options);
[dP_2 R_2 C_2 negLL exitflag] = PAL_SDT_ROCML_Fit(cumNumHF_2, OutOfNum_2,SDTF,invSDTF,'searchOptions',options);
[dP_3 R_3 C_3 negLL exitflag] = PAL_SDT_ROCML_Fit(cumNumHF_3, OutOfNum_3,SDTF,invSDTF,'searchOptions',options);

message=sprintf('\rParameters d-prime, ratio-of-SDs and four criteria C now fitted for each ROC curve:');
disp(message);

% Calculate Z scores from cumulative proportions
ZHF_1 = PAL_PtoZ(pHF_1);  
ZHF_2 = PAL_PtoZ(pHF_2); 
ZHF_3 = PAL_PtoZ(pHF_3); 

%Create continuous ROC curves based on the fitted parameter estimates
vecC = [-2.5:0.1:2.5];

vecDP_1 = ones(1,length(vecC)).*dP_1;
vecDP_2 = ones(1,length(vecC)).*dP_2;
vecDP_3 = ones(1,length(vecC)).*dP_3;

vecR_1 = ones(1,length(vecC)).*R_1;
vecR_2 = ones(1,length(vecC)).*R_2;
vecR_3 = ones(1,length(vecC)).*R_3;

linePHF_1 = invSDTF(vecDP_1,vecC,'ratioSDvalue',vecR_1);
linePHF_2 = invSDTF(vecDP_2,vecC,'ratioSDvalue',vecR_2);
linePHF_3 = invSDTF(vecDP_3,vecC,'ratioSDvalue',vecR_3);

lineZHF_1 = PAL_PtoZ(linePHF_1);
lineZHF_2 = PAL_PtoZ(linePHF_2);
lineZHF_3 = PAL_PtoZ(linePHF_3);

% Create 1st figure for curves of prop. hits versus false alarms
figure('name','Prop. Hits v. False alarms, and z(H) v. z(F), for various d-prime and SD ratios','units','pixels',...
    'position',[50 50 1000 400]);
subplot(1,2,1);
set(gca, 'fontsize',16);
a = plot(linePHF_1(:,2),linePHF_1(:,1),'g-','linewidth',2);
hold on;
b = plot(pHF_1(:,2),pHF_1(:,1), 'o','MarkerEdgeColor','k','MarkerFaceColor','k','Markersize',8);
hold on;
plot(linePHF_2(:,2),linePHF_2(:,1),'g-','linewidth',2), axis square;
hold on;
c = plot(pHF_2(:,2),pHF_2(:,1), 's','MarkerEdgeColor','k','MarkerFaceColor','k','Markersize',8);
hold on;
plot(linePHF_3(:,2),linePHF_3(:,1),'g-','linewidth',2), axis square;
hold on;
d = plot(pHF_3(:,2),pHF_3(:,1), '*','MarkerEdgeColor','k','MarkerFaceColor','k','Markersize',8);
hold on;

% Set graph properties
xlim([0.0,1.0]);
ylim([0.0,1.0]);
set(gca, 'Xtick',[0.0:0.2:1.0]);
set(gca, 'Ytick',[0.0:0.2:1.0]);
ylabel('Proportion hits','FontSize',16);
xlabel('Proportion false alarms','FontSize',16);

% Add in legend
DPstrng1 = num2str(dP_1,'%5.2f');
DPstrng2 = num2str(dP_2,'%5.2f');
DPstrng3 = num2str(dP_3,'%5.2f');
Rstrng1 = num2str(R_1,'%5.2f');
Rstrng2 = num2str(R_2,'%5.2f');
Rstrng3 = num2str(R_3,'%5.2f');
string1 = strcat('dp=',DPstrng1,'  R=',Rstrng1);
string2 = strcat('dp=',DPstrng2,'  R=',Rstrng2);
string3 = strcat('dp=',DPstrng3,'  R=',Rstrng3);
h = legend([a b c d],'fits',string1,string2,string3,2);
set(h,'Interpreter','none','fontsize',16,'Location','SouthEast');

% Create 2nd figure for curves of z(H) versus z(F)
subplot(1,2,2);
set(gca, 'fontsize',16);
a = plot(lineZHF_1(:,2),lineZHF_1(:,1),'g-','linewidth',2);
hold on;
b = plot(ZHF_1(:,2),ZHF_1(:,1), 'o','MarkerEdgeColor','k','MarkerFaceColor','k','Markersize',8);
hold on;
plot(lineZHF_2(:,2),lineZHF_2(:,1),'g-','linewidth',2), axis square;
hold on;
c = plot(ZHF_2(:,2),ZHF_2(:,1), 's','MarkerEdgeColor','k','MarkerFaceColor','k','Markersize',8);
hold on;
plot(lineZHF_3(:,2),lineZHF_3(:,1),'g-','linewidth',2), axis square;
hold on;
d = plot(ZHF_3(:,2),ZHF_3(:,1), '*','MarkerEdgeColor','k','MarkerFaceColor','k','Markersize',8);

% Set graph properties
xlim([-2.5,2.5]);
ylim([-2.5,2.5]);
set(gca, 'Xtick',[-2.5:.5:2.5]);
set(gca, 'Ytick',[-2.5:.5:2.5]);
ylabel('z(H)','FontSize',16);
xlabel('z(F)','FontSize',16);
drawnow

message = sprintf('\r----------First ROC---------');
disp(message);

disp('Determining standard errors.....');
paramsValues = [dP_1 R_1 C_1];
[SD_dP_1 SD_R_1 paramsSim LLsim converged] = PAL_SDT_ROCML_BootstrapParametric(paramsValues, OutOfNum_1, SDTF, invSDTF, Bse,'searchOptions',options);

disp('Determining goodness-of-fit.....');
[Dev pDev] = PAL_SDT_ROCML_GoodnessOfFit(paramsValues, cumNumHF_1, OutOfNum_1, SDTF, invSDTF, Bgof,'searchOptions',options);

disp('Determining if ratio-of-SDs different from 1 .....');
[TLR pTLR] = PAL_SDT_ROCML_RatioSDcomparison(cumNumHF_1, OutOfNum_1, 1, SDTF, invSDTF, Brsd,'searchOptions',options);

message = sprintf('\rFull results:');
disp(message);
message = sprintf('d-prime: %6.4f',dP_1);
disp(message);
message = sprintf('SE d-prime: %6.4f',SD_dP_1);
disp(message);
message = sprintf('ratioSD: %6.4f',R_1);
disp(message);
message = sprintf('SE ratioSD: %6.4f',SD_R_1);
disp(message);
message = sprintf('Deviance: %6.4f',Dev);
disp(message);
message = sprintf('Goodness-of-fit p value: %6.4f',pDev);
disp(message);
message = sprintf('Ratio-of-SDs-different-from-1 p value: %6.4f\r',pTLR);
disp(message);


disp('--------Second ROC----------');

disp('Determining standard errors.....');
paramsValues = [dP_2 R_2 C_2];
[SD_dP_2 SD_R_2 paramsSim LLsim converged] = PAL_SDT_ROCML_BootstrapParametric(paramsValues, OutOfNum_2, SDTF, invSDTF, Bse,'searchOptions',options);

disp('Determining goodness-of-fit....');
[Dev pDev] = PAL_SDT_ROCML_GoodnessOfFit(paramsValues, cumNumHF_2, OutOfNum_2, SDTF, invSDTF, Bgof,'searchOptions',options);

disp('Determining if ratio-of-SDs different from 1 ....');
[TLR pTLR] = PAL_SDT_ROCML_RatioSDcomparison(cumNumHF_2, OutOfNum_2, 1, SDTF, invSDTF, Brsd,'searchOptions',options);

message = sprintf('\rFull results:');
disp(message);
message = sprintf('d-prime: %6.4f',dP_2);
disp(message);
message = sprintf('SE d-prime: %6.4f',SD_dP_2);
disp(message);
message = sprintf('ratioSD: %6.4f',R_2);
disp(message);
message = sprintf('SE ratioSD: %6.4f',SD_R_2);
disp(message);
message = sprintf('Deviance: %6.4f',Dev);
disp(message);
message = sprintf('Goodness-of-fit p value: %6.4f',pDev);
disp(message);
message = sprintf('Ratio-of-SDs-different-from-1 p value: %6.4f\r',pTLR);
disp(message);


disp('---------Third ROC------------');

disp('Determining standard errors....');
paramsValues = [dP_3 R_3 C_3];
[SD_dP_3 SD_R_3 paramsSim LLsim converged] = PAL_SDT_ROCML_BootstrapParametric(paramsValues, OutOfNum_3, SDTF, invSDTF, Bse,'searchOptions',options);

disp('Determining goodness-of-fit....');
[Dev pDev] = PAL_SDT_ROCML_GoodnessOfFit(paramsValues, cumNumHF_3, OutOfNum_3, SDTF, invSDTF, Bgof,'searchOptions',options);

disp('Determining if ratio-of-SDs different from 1 ....');
[TLR pTLR] = PAL_SDT_ROCML_RatioSDcomparison(cumNumHF_3, OutOfNum_3, 1, SDTF, invSDTF, Brsd,'searchOptions',options);

message=sprintf('\rFull results:');
disp(message);
message = sprintf('d-prime: %6.4f',dP_3);
disp(message);
message = sprintf('SE d-prime: %6.4f',SD_dP_3);
disp(message);
message = sprintf('ratio SD: %6.4f',R_3);
disp(message);
message = sprintf('SE ratioSD: %6.4f',SD_R_3);
disp(message);
message = sprintf('Deviance: %6.4f',Dev);
disp(message);
message = sprintf('Goodness-of-fit p value: %6.4f',pDev);
disp(message);
message = sprintf('Ratio-of-SDs-different-from-1 p value: %6.4f\r',pTLR);
disp(message);
