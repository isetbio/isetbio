%
% PAL_SDT_PSvAS_SummSquare_Demo
% Demonstrates use of Palamedes functions to determine whether the 
% data from a 5-PF (psychometric function) summation square 
% experiment, for the detection of two stimuli in the target interval, 
% accords more with probability summation (PS) or additive summation (AS) 
% under the assumptions of signal-detection-theory (SDT) and assuming
% that the observer is monitoring both channels sensitive to the two
% stimuli
%
% Note that the Bootstrap and Goodness-of-fit routines in this script will
% likely take some time to execute depending on the number of specified 
% simulations and the speed of the computer.  An indication of the time 
% to execute on a MacBook Pro running under OSX 10.9.5 is given when 
% the script is executed
%
%Demonstrates usage of Palamedes SDT summation multiple psychometric 
% function (PF) fitting routines:
%-PAL_SDT_Summ_MultiplePFML_Fit
%-PAL_SDT_Summ_MultiplePFML_BootstrapParametric
%-PAL_SDT_Summ_MultiplePFML_GoodnessOfFit
%
%Demonstrates usage of Palamedes SDT PS (probability) and AS (additive)
% summation routines:
%-PAL_SDT_PS_uneqSLtoPC
%-PAL_SDT_AS_uneqSLtoPC
%
%Demonstrates usage of Palamedes PF fitting routines:
%-PAL_PFML_Fit;
%
%Demonstrates usage of PF routines, e.g.
%-PAL_Logistic
%
%More information on any of these functions may be found by typing
% help followed by the name of the function, 
% e.g., help PAL_SDT_AS_uneqSLtoPC
%
%Introduced: Palamedes version 1.8.0 (FK & NP)

clear all;

message = sprintf('Number of simulations to determine standard errors: ');
Bse = input(message);
message = sprintf('Number of simulations to determine Goodness-of-Fit p-values: ');
Bmc = input(message);

MacProTimePerSim=0.5;
timeToExecute=round((Bse+Bmc)*MacProTimePerSim);
message=sprintf('\rApprox. time to execute simulations on a MacBook Pro OSX 10.9.5 is %4d mins',timeToExecute); 
disp(message);


% Stimulus Levels for 5 PFs, include zeros for 'empty' monitored channels
StimLevels(1,1,:)=[1 2 3 4 5]; % stim A 
StimLevels(1,2,:)=[0 0 0 0 0]; % stim B for channel B 

StimLevels(2,1,:)=[0.8 1.6 2.4 3.2 4.0]; % stim A 
StimLevels(2,2,:)=[0.2 0.4 0.6 0.8 1]; % stim B 

StimLevels(3,1,:)=[0.5 1 1.5 2 2.5]; % stim A 
StimLevels(3,2,:)=[0.5 1 1.5 2 2.5]; % stim B 

StimLevels(4,1,:)=[0.2 0.4 0.6 0.8 1]; % stim A 
StimLevels(4,2,:)=[0.8 1.6 2.4 3.2 4.0]; % stim B

StimLevels(5,1,:)=[0 0 0 0 0]; % stim A for channel B
StimLevels(5,2,:)=[1 2 3 4 5]; % for stim B 


% Number of correct responses for each stimulus level pair in the 5 PFs
NumPos(1,:)=[50 63 77 86 94];
NumPos(2,:)=[52 60 70 83 94]; 
NumPos(3,:)=[50 58 68 80 92];
NumPos(4,:)=[51 62 74 85 96];
NumPos(5,:)=[54 69 82 92 99];


%Number of trials for each stimulus level pair in the 5 PFs 
OutOfNum(1,:) = [100 100 100 100 100];
OutOfNum(2,:) = [100 100 100 100 100];
OutOfNum(3,:) = [100 100 100 100 100];
OutOfNum(4,:) = [100 100 100 100 100];
OutOfNum(5,:) = [100 100 100 100 100];



%------- Fit each PF with a standard model, such as a Logistic or Wibull,
% in order to estimate thresholds that can be put onto the graph.  This is 
% for illustrative purposes only - the fitted parameters here are NOT used 
% in the modeling of probability and additiv summation

searchGrid.alpha = 0:.001:5; % range of possible threshold values
searchGrid.beta = logspace(.1,3,100); % range of possible slope values
searchGrid.gamma = .5;  %guessing rate
searchGrid.lambda = 0.0;  %lapse rate

%Threshold and Slope are free parameters, guess and lapse rate are fixed
paramsFree = [1 1 0 0];  %1: free parameter, 0: fixed parameter
 
%Fit Logistic functions
PF = @PAL_Logistic;  %Alternatives: PAL_Gumbel, PAL_Weibull, 
                     %PAL_CumulativeNormal, PAL_HyperbolicSecant

[nrows, ncols, nnums]=size(StimLevels);

% Fit each data PF with a standard model in order to obtain thresholds
% to show on the graph
threshX=zeros(1,nrows);
threshY=zeros(1,nrows);

for i=1:nrows
    if (mean(StimLevels(i,1,:))==0)
        threshX(i)=0;
    else
    Stim(1,:)=StimLevels(i,1,:);
    [params, trash, trash, trash] = PAL_PFML_Fit(Stim(1,:),NumPos(i,:),OutOfNum(i,:),searchGrid,paramsFree,PF,'searchOptions',[]);
    threshX(i)=params(1);
    end
    
    if (mean(StimLevels(i,2,:))==0)
        threshY(i)=0;
    else
    Stim(1,:)=StimLevels(i,2,:);
    [params, trash, trash, trash] = PAL_PFML_Fit(Stim(1,:),NumPos(i,:),OutOfNum(i,:),searchGrid,paramsFree,PF,'searchOptions',[]);
    threshY(i)=params(1);
    end
end

% Plot data thresholds
figure('name','Individual Psychometric Functions','units','pixels',...
    'position',[50 50 400 400]); 

a=plot(threshX,threshY,'ko','markersize',8,'markerfacecolor','k');
hold on;
set(gca, 'fontsize',18);
set(gca, 'Xtick',0:0.5:3.5);
set(gca, 'Ytick',0:0.5:3.5);
axis([0 3.5 0 3.5]); axis square;
xlabel('Stimulus A');
ylabel('Stimulus B');
title('Summation square');
drawnow


%--------- Fit PFs with probability summation (PS) model ----------

% Fit PS model to stim A, B and the combinations of A+B in the 5 PFs, 
% simultaneously, in order to obtain parameters gA, gB, pA and pB, 
% where g is stimulus gain and p the transducer exponent.  The fitting 
% procedure also provides the negative log likeleihood of the fit which
% is usd to compare the PS with the AS model

Q=2; % Number of monitored channels
M=2; % Number of alternatives in forced-choice task

SummFunc=@PAL_SDT_PS_uneqSLtoPC; % use PS model with unequal stimulus 
% levels to convert stimulus intensity to proportion correct

% Initial guesses for gA, gB, pA and pB
PSgParams=[0.3 0.3];
PSpParams=[1.5 1.5];

% Fit all five data PFs simultaneously with PS model
[PSgParams, PSpParams, PSnegLL, trash, trash] = PAL_SDT_Summ_MultiplePFML_Fit(StimLevels,PSgParams,PSpParams,NumPos,OutOfNum,SummFunc,M,Q);


%----------- Fit PFs with additive summation (AS) model ------------------

% Same as above for addiitve summation model

SummFunc=@PAL_SDT_AS_uneqSLtoPC;

% Initial guesses for gA, gB, pA and pB
ASgParams=[0.3 0.3];
ASpParams=[1.5 1.5];

% Fit all PF data simultaneously with AS summation model
[ASgParams, ASpParams, ASnegLL, exitflag, output] = PAL_SDT_Summ_MultiplePFML_Fit(StimLevels,ASgParams,ASpParams,NumPos,OutOfNum,SummFunc,M,Q);



% ------ Generate new PFs from fitted summation model and fit these with the same model as for the original data
% This is for illustrative purposes in the graph and is not part of the modelling itself -----

% ---- First for PS model---
StimLevelsBaseline=StimLevels(1,1,:);
meanOutNum=mean(mean(OutOfNum));
newOutNum=ones(1,nnums).*meanOutNum;
newPC=zeros(1,nnums);

numModelPoints=100;
minLogR=-2.5;
maxLogR=2.5;
incLogR=(maxLogR-minLogR)/numModelPoints;

PSthresh=zeros(2,numModelPoints);

for k=1:numModelPoints
    r=minLogR+k*incLogR;
    r=10.^r;
    x1=(r.*StimLevelsBaseline)./(1+r);
    x2=StimLevelsBaseline-x1;
   
    for i=1:nnums
        newPC(i)=PAL_SDT_PS_uneqSLtoPC([x1(i) x2(i)],PSgParams,PSpParams,M,Q);
    end
    
    Pos=newPC.*meanOutNum;
    
    [params, trash, trash, trash] = PAL_PFML_Fit(x1,Pos,newOutNum,searchGrid,paramsFree,PF,'searchOptions',[]);
    PSthresh(1,k)=params(1);
    [params, trash, trash, trash] = PAL_PFML_Fit(x2,Pos,newOutNum,searchGrid,paramsFree,PF,'searchOptions',[]);
    PSthresh(2,k)=params(1);
    
end


% -----Second for the AS model -----

ASthresh=zeros(2,numModelPoints);

for k=1:numModelPoints
    r=minLogR+k*incLogR;
    r=10.^r;
    x1=(r.*StimLevelsBaseline)./(1+r);
    x2=StimLevelsBaseline-x1;
    
    for i=1:nnums
        newPC(i)=PAL_SDT_AS_uneqSLtoPC([x1(i) x2(i)],ASgParams,ASpParams,M,Q);
    end
    
    Pos=newPC.*meanOutNum;
    
    [params, trash, trash, trash] = PAL_PFML_Fit(x1,Pos,newOutNum,searchGrid,paramsFree,PF,'searchOptions',[]);
    ASthresh(1,k)=params(1);
    [params, trash, trash, trash] = PAL_PFML_Fit(x2,Pos,newOutNum,searchGrid,paramsFree,PF,'searchOptions',[]);
    ASthresh(2,k)=params(1);
   
end


% Plot model predictions on graph
b=plot(PSthresh(1,:),PSthresh(2,:),'-','linewidth',2,'color',[0 .9 0]);
hold on;
c=plot(ASthresh(1,:),ASthresh(2,:),'-','linewidth',2,'color',[.9 0 0]);

% Add in legend
h = legend([a b c],'Data','PS Model','AS model');
set(h,'Interpreter','none','fontsize',16,'Location','NorthEast');
drawnow




%---------Determine bootstrap errors on fitted g and p parameters--------



SummFunc=@PAL_SDT_PS_uneqSLtoPC;
fprintf('\rConducting bootstrap simulations for Prob. Summ. (PS) model:   ')
[PSgSE, PSpSE] = PAL_SDT_Summ_MultiplePFML_BootstrapParametric(StimLevels,PSgParams,PSpParams,OutOfNum,SummFunc,M,Q,Bse);

SummFunc=@PAL_SDT_PS_uneqSLtoPC;
fprintf('\nConducting bootstrap simulations for Add. Summ. (AS) model:   ')
[ASgSE, ASpSE] = PAL_SDT_Summ_MultiplePFML_BootstrapParametric(StimLevels,ASgParams,ASpParams,OutOfNum,SummFunc,M,Q,Bse);


% ------------ Type out parameter estimates for both models-----------

message=sprintf('\rEstimated PS model params for stimulus A: gA: %6.4f se: %6.4f,  pA: %6.4f se: %6.4f',PSgParams(1),PSgSE(1),PSpParams(1),PSpSE(1)); 
disp(message);
message=sprintf('Estimated PS model params for stimulus B: gB: %6.4f se: %6.4f,  pB: %6.4f se: %6.4f',PSgParams(2),PSgSE(2),PSpParams(2),PSpSE(2)); 
disp(message);
message=sprintf('Estimated AS model params for stimulus A: gA: %6.4f sd: %6.4f,  pA: %6.4f se: %6.4f',ASgParams(1),ASgSE(1),ASpParams(1),ASpSE(1)); 
disp(message);
message=sprintf('Estimated AS model params for stimulus B: gB: %6.4f sd: %6.4f,  pB: %6.4f se: %6.4f',ASgParams(2),ASgSE(2),ASpParams(2),ASpSE(2)); 
disp(message);



%---------------------- Model comparisons ------------------------------

% AS-PS difference between Akaike's AIC.  Note that a positive value 
% means that PS model is better, a negative value that AS model is better
PSaic = -2*(-PSnegLL)+2*4; % Note four free parameters per model
ASaic = -2*(-ASnegLL)+2*4;
deltaAIC = ASaic - PSaic;
message = sprintf('\rDifference between models in Akaike Information Criterion (AIC) is: %8.4f',deltaAIC);
disp(message);
message = sprintf('Note that positive value = PS better, negative = AS better');
disp(message);


%-----------------------Goodness of fits-------------------------------

SummFunc=@PAL_SDT_PS_uneqSLtoPC;
fprintf('Conducting Goodness-of-fit (GOF) simulations for PS model:   ')
[trash, pDev, trash, trash] = PAL_SDT_Summ_MultiplePFML_GoodnessOfFit(StimLevels,PSgParams,PSpParams,NumPos,OutOfNum,SummFunc,M,Q,Bmc);
PSpDev=pDev;

SummFunc=@PAL_SDT_AS_uneqSLtoPC;
fprintf('\nConducting Goodness-of-fit (GOF) simulations for AS model:   ')
[trash, pDev, DevSim, converged] = PAL_SDT_Summ_MultiplePFML_GoodnessOfFit(StimLevels,ASgParams,ASpParams,NumPos,OutOfNum,SummFunc,M,Q,Bmc);
ASpDev=pDev;

message=sprintf('\rGOF p-values (e.g. p<0.05 to reject model): PS p=%6.4f;  AS p=%6.4f',PSpDev,ASpDev); 
disp(message);

