%
% PAL_SDT_PSvAS_3PFmultipleFit_Demo
% Demonstrates use of Palamedes functions to determine whether the 
% psychometric functions for the detection of two stimuli in the target 
% interval accords with probability or additive summation 
%
% Note that the Bootstrap and Goodness-of-fit routines in this script will
% likely take some time to execute depending on the number of specified 
% simulations and the speed of the computer.  An indication of the time 
% to execute on a MacBook Pro running under OSX 10.9.5 is given when 
% the script is executed
%
%Demonstrates usage of Palamedes SDT summation PF fitting routines:
%-PAL_SDT_Summ_MultiplePFML_Fit
%-PAL_SDT_Summ_MultiplePFML_BootstrapParametric
%-PAL_SDT_Summ_MultiplePFML_GoodnessOfFit
%
%Demonstrates usage of Palamedes SDT probability summation (PS) routines:
%-PAL_SDT_PS_PCtoSL
%-PAL_SDT_PS_uneqSLtoPC
%
%Demonstrates usage of Palamedes SDT additive summation (AS) routines:
%-PAL_SDT_AS_PCtoSL
%-PAL_SDT_AS_uneqSLtoPC
%
%
%More information on any of these functions may be found by typing
%help followed by the name of the function. e.g., help PAL_SDT_AS_PCtoSL
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

Q=2; % Number of monitored channels
M=2; % Number of alternatives in forced-choice task

% Stimulus Levels for 3PFs, include zeros for empty monitored channels
StimLevels(1,1,:)=[1 2 3 4 5 6 7 8]; % for stim A alone
StimLevels(1,2,:)=[0 0 0 0 0 0 0 0]; % channel B for stim A alone
StimLevels(2,1,:)=[0 0 0 0 0 0 0 0]; % channel A for stim B alone
StimLevels(2,2,:)=[1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5]; %Stim B alone
StimLevels(3,1,:)=[1.2 2.2 3.2 4.2 5.2 6.2 7.2 8.2]; %Stim A in Stim A+B
StimLevels(3,2,:)=[1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5]; %Stim B in Stim A+B

%for displaying Stim A+B data
StimLevelsABindex(1,:) = (StimLevels(3,1,:) + StimLevels(3,2,:))./2; 

% Number correct reponses for A, B and A+B stimuli
NumPos(1,:)=[44 49 53 60 64 72 73 78];
NumPos(2,:)=[43 46 56 58 68 70 74 79]; 
NumPos(3,:)=[44 57 59 68 76 76 79 79];

%Number of trials for A, B and A+B. May contain zeros. 
OutOfNum(1,:) = [80 80 80 80 80 80 80 80];
OutOfNum(2,:) = [80 80 80 80 80 80 80 80];
OutOfNum(3,:) = [80 80 80 80 80 80 80 80];


%Plot raw data, proportion correct PC against stimulus level
PC = NumPos ./ OutOfNum; 

figure('name','Individual Psychometric Functions','units','pixels',...
    'position',[50 50 1000 400]); 


subplot(1,2,1);
axis([1 8 .4 1.0]); axis square;
set(gca, 'Xtick',1:1:10);
set(gca, 'Ytick',.4:.1:1);
m(1,1,:)=StimLevels(1,1,:);
StimLevelsPlot=squeeze(m);
a=plot(StimLevelsPlot,PC(1,:),'ro','markersize',8,'markerfacecolor','r'); 
hold on;
m(1,1,:)=StimLevels(2,2,:);
StimLevelsPlot=squeeze(m);
b=plot(StimLevelsPlot,PC(2,:),'bo','markersize',8,'markerfacecolor','b');
hold on;
c=plot(StimLevelsABindex,PC(3,:),'go','markersize',8,'markerfacecolor','g');
hold on;
set(gca, 'fontsize',18);
xlabel('Stimulus level');
ylabel('Proportion correct');
title('Probability summation');

subplot(1,2,2);

axis([1 8 .4 1.0]); axis square;
set(gca, 'Xtick',1:1:10);
set(gca, 'Ytick',.4:.1:1);
m(1,1,:)=StimLevels(1,1,:);
StimLevelsPlot=squeeze(m);
plot(StimLevelsPlot,PC(1,:),'ro','markersize',8,'markerfacecolor','r'); 
hold on;
m(1,1,:)=StimLevels(2,2,:);
StimLevelsPlot=squeeze(m);
plot(StimLevelsPlot,PC(2,:),'bo','markersize',8,'markerfacecolor','b');
hold on;
plot(StimLevelsABindex,PC(3,:),'go','markersize',8,'markerfacecolor','g');
hold on;
set(gca, 'fontsize',18);
xlabel('Stimulus level');
ylabel('Proportion correct');
title('Additive summation');
drawnow


%------Calculate PF prediction for probability summation (PS) model-------

% Fit PS model to stim A, B and A+B simultaneously to obtain parameters 
% gA, gB, pA and pB, where g is stimulus gain and p transducer exponent
SummFunc=@PAL_SDT_PS_uneqSLtoPC;

% Initial guesses for gA, gB, pA and pB
PSgParams=[0.3 0.3];
PSpParams=[1.5 1.5];

% Fit all three data functions simultaneously with PS model
[PSgParams, PSpParams, PSnegLL, exitflag, output] = PAL_SDT_Summ_MultiplePFML_Fit(StimLevels,PSgParams,PSpParams,NumPos,OutOfNum,SummFunc,M,Q);


% Now create continuous functions to show model fits to data
StimLevelsFineGrain(1,:) = linspace(min(StimLevels(1,1,:)),max(StimLevels(1,1,:)),1000);
StimLevelsFineGrain(2,:) = linspace(min(StimLevels(2,2,:)),max(StimLevels(2,2,:)),1000);
StimLevelsFineGrain(3,:) = linspace(min(StimLevels(3,1,:)),max(StimLevels(3,1,:)),1000);
StimLevelsFineGrain(4,:) = linspace(min(StimLevels(3,2,:)),max(StimLevels(3,2,:)),1000);
StimLevelsFineGrainABindex = linspace(min(StimLevelsABindex),max(StimLevelsABindex),1000);

PCfromPSfineGrain(1,:) = PAL_SDT_PS_SLtoPC(StimLevelsFineGrain(1,:),PSgParams(1),PSpParams(1),M,Q,1);
PCfromPSfineGrain(2,:) = PAL_SDT_PS_SLtoPC(StimLevelsFineGrain(2,:),PSgParams(2),PSpParams(2),M,Q,1);

for i=1:length(StimLevelsFineGrainABindex)
PCfromPSfineGrain(3,i) = PAL_SDT_PS_uneqSLtoPC([StimLevelsFineGrain(3,i) StimLevelsFineGrain(4,i)],PSgParams,PSpParams,M,Q);
end

subplot(1,2,1);
d=plot(StimLevelsFineGrain(1,:),PCfromPSfineGrain(1,:),'r-','linewidth',2);
hold on;
e=plot(StimLevelsFineGrain(2,:),PCfromPSfineGrain(2,:),'b-','linewidth',2);
hold on;
f=plot(StimLevelsFineGrainABindex,PCfromPSfineGrain(3,:),'g-','linewidth',2);
drawnow
hold on;

h = legend([a b c d e f],'Data A','Data B','Data A+B','Model A','Model B','Model A+B');
set(h,'Interpreter','none','fontsize',16,'Location','SouthEast');
drawnow

%-----------Calculate PF for additive summation model----------------------

% Fit AS model to stim A, B and A+B simultaneously to obtain parameters gA, gB, pA and pB
SummFunc=@PAL_SDT_AS_uneqSLtoPC;

% Initial guesses for gA, gB, pA and pB
ASgParams=[0.3 0.3];
ASpParams=[1.5 1.5];

% Fit all three data functions simultaneously with AS model
[ASgParams, ASpParams, ASnegLL, trash, trash] = PAL_SDT_Summ_MultiplePFML_Fit(StimLevels,ASgParams,ASpParams,NumPos,OutOfNum,SummFunc,M,Q);


% Now create continuous functions to show model fits to data
PCfromASfineGrain(1,:) = PAL_SDT_AS_SLtoPC(StimLevelsFineGrain(1,:),ASgParams(1),ASpParams(1),M,Q,1);
PCfromASfineGrain(2,:) = PAL_SDT_AS_SLtoPC(StimLevelsFineGrain(2,:),ASgParams(2),ASpParams(2),M,Q,1);

for i=1:length(StimLevelsFineGrainABindex)
PCfromASfineGrain(3,i) = PAL_SDT_AS_uneqSLtoPC([StimLevelsFineGrain(3,i) StimLevelsFineGrain(4,i)],ASgParams,ASpParams,M,Q);
end

subplot(1,2,2);
d=plot(StimLevelsFineGrain(1,:),PCfromASfineGrain(1,:),'r-','linewidth',2);
hold on;
e=plot(StimLevelsFineGrain(2,:),PCfromASfineGrain(2,:),'b-','linewidth',2);
hold on;
f=plot(StimLevelsFineGrainABindex,PCfromASfineGrain(3,:),'g-','linewidth',2);
drawnow
hold on;

%---------Determine bootstrap errors on fitted g and p parameters--------

SummFunc=@PAL_SDT_PS_uneqSLtoPC;
fprintf('\rConducting simulations for Prob. Summ. (PS) model:   ')
[PSgSE, PSpSE] = PAL_SDT_Summ_MultiplePFML_BootstrapParametric(StimLevels,PSgParams,PSpParams,OutOfNum,SummFunc,M,Q,Bse);

SummFunc=@PAL_SDT_PS_uneqSLtoPC;
fprintf('\nConducting simulations for Add. Summ. (AS) model:   ')
[ASgSE, ASpSE] = PAL_SDT_Summ_MultiplePFML_BootstrapParametric(StimLevels,ASgParams,ASpParams,OutOfNum,SummFunc,M,Q,Bse);


% ------------ Type out parameter estimates for both models-----------

message=sprintf('\rEstimated PS model params for stimulus A: gA: %6.4f se: %6.4f,  pA: %6.4f se: %6.4f',PSgParams(1),PSgSE(1),PSpParams(1),PSpSE(1)); 
disp(message);
message=sprintf('Estimated PS model params for stimulus B: gB: %6.4f se: %6.4f,  pB: %6.4f se: %6.4f',PSgParams(2),PSgSE(2),PSpParams(2),PSpSE(2)); 
disp(message);
message=sprintf('Estimated AS model params for stimulus A: gA: %6.4f se: %6.4f,  pA: %6.4f se: %6.4f',ASgParams(1),ASgSE(1),ASpParams(1),ASpSE(1)); 
disp(message);
message=sprintf('Estimated AS model params for stimulus B: gB: %6.4f se: %6.4f,  pB: %6.4f se: %6.4f',ASgParams(2),ASgSE(2),ASpParams(2),ASpSE(2)); 
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


%----- Goodness of fits-----------

SummFunc=@PAL_SDT_PS_uneqSLtoPC;
fprintf('Conducting Goodness-of-fit (GOF) simulations for PS model:   ')
[trash, pDev, trash, trash] = PAL_SDT_Summ_MultiplePFML_GoodnessOfFit(StimLevels,PSgParams,PSpParams,NumPos,OutOfNum,SummFunc,M,Q,Bmc);
PSpDev=pDev;

SummFunc=@PAL_SDT_AS_uneqSLtoPC;
fprintf('\nConducting Goodness-of-fit (GOF) simulations for AS model:   ')
[trash, pDev, trash, trash] = PAL_SDT_Summ_MultiplePFML_GoodnessOfFit(StimLevels,ASgParams,ASpParams,NumPos,OutOfNum,SummFunc,M,Q,Bmc);
ASpDev=pDev;

message=sprintf('\rGOF p-values (e.g. p<0.05 to reject model): PS p=%6.4f;  AS p=%6.4f',PSpDev,ASpDev); 
disp(message);

