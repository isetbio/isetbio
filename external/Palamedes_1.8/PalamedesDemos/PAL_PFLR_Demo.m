%
%PAL_PFLR_Demo  Demonstrates use of Palamedes functions to (1) fit
%psychometric functions to multiple conditions simultaneously (2) determine
%standard errors of the free parameters in simultaneous fits using
%bootstrap simulations and (3) perform statistical model comparisons 
%between nested models. Different model comparisons are performed including
%a goodness-of-fit test for a user-defined model.
%
%Program first fits psychometric functions (PFs) to some data in a two
%group experiment. The model fits individual thresholds and slopes to the 
%two conditions, assumes that the guess rate is 0.5 and the lapse rate 
%equals 0. Thus, this is a 4 parameter model (2 thresholds, 2 slopes). 
%Plots of the fit are presented in a figure.
%
%Program then estimates standard errors on the 4 free parameters using
%bootstrap simulations. Standard error bars are added to the figure.
%
%Program then performs 4 model comparisons. Figure 2 shows results of a
%model comparison between the above model and a model which assumes that
%thresholds and slopes are equal between conditions. Figure 3 shows results
%of a comparison between the above model and one that assumes that
%thresholds are equal but allows slopes to vary. Figure 4 compares above
%model to model which assumes slopes are equal but allows thresholds to
%vary. Figure 5 compares a model which assumes that slopes are equal but
%allows thresholds to vary to the saturated model (i.e., this is a
%Goodness-of-Fit test).
%
%Histograms of the simulated transformed likelihood ratios 
%[TLR = -2x(log likelihood lesser model - log likelihood fuller model)] are 
%created. Also displayed for each model comparison are the proportion of 
%simulated TLR values that are greater than the data's TLR value.
%
%In case the Matlab routines chi2cdf and chi2pdf are present, program
%superimposes theoretical chi-square distribution with appropriate degrees 
%of freedom on empirical sampling distributions and presents the p-value 
%based on the chi-square distribution with the appropriate degrees of 
%freedom.
%
%Demonstrates usage of Palamedes functions:
%-PAL_PFML_FitMultiple
%-PAL_PFML_BootstrapParametricMultiple
%-PAL_PFML_BootstrapNonParametricMultiple
%-PAL_PFML_GoodnessOfFitMultiple
%-PAL_PFLR_ModelComparison
%secondary:
%PAL_Logistic
%
%More information on any of these functions may be found by typing
%help followed by the name of the function. e.g., help PAL_PFML_FitMultiple
%
%NP (September 2009)
%NP Modified to utilize new output argument of PAL_PFML_FitMultiple: 
%   numParams (November 2009)

clear all;

if exist('OCTAVE_VERSION');
    fprintf('\nUnder Octave, Figures do not render (at all!) as intended. Visit\n');
    fprintf('www.palamedestoolbox.org/demosfiguregallery.html to see figures\n');
    fprintf('as intended.\n\n');
    cont = input('Do you wish to continue? [y/n]: ', 's');
    if strcmp(cont,'y')
        fprintf('\n\n');
    else
        break;
    end
end

message = 'Parametric Bootstrap (1) or Non-Parametric Bootstrap? (2): ';
ParOrNonPar = input(message);
message = sprintf('Number of simulations to perform to determine standar');
message = strcat(message, 'd errors: ');
Bse = input(message);
message = sprintf('Number of simulations to perform to determine model c');
message = strcat(message, 'omparison p-values: ');
Bmc = input(message);

tic

%Stimulus Levels 2 (conditions) x 5 (stimulus levels)
%rows need not be identical
StimLevels = [-2:1:2; -2:1:2];

%Number of trials at each entry of StimLevels
%entries may be zero (useful if stimulus levels are not identical between
%conditions).
OutOfNum = [100 100 100 100 100; 100 100 100 100 100];

%Number of positive responses (e.g., 'yes', 'correct') for each entry of 
%StimLevels
NumPos = [61 70 81 92 97; 59 59 67 86 91];

%Plot raw data
ProportionCorrectObserved = NumPos ./ OutOfNum; 
StimLevelsFineGrain = [min(min(StimLevels)):(max(max(StimLevels) - ... 
    min(min(StimLevels))))./1000:max(max(StimLevels))];
figure('name','Individual Psychometric Functions','units','pixels',...
    'position',[100 100 500 500]);
plot(StimLevels(1,:),ProportionCorrectObserved(1,:),'ko','markersize',...
    10,'markerfacecolor','k');
h1 = gca;
set(h1, 'units','pixels','position',[75 300 375 175]);
set(h1, 'fontsize',16);
set(h1, 'Xtick',StimLevels(1,:));
set(h1, 'Ytick',[.5:.1:1]);
axis([-2.5 2.5 .5 1]);
hold on;
plot(StimLevels(2,:),ProportionCorrectObserved(2,:),'ks','markersize',...
    10,'markerfacecolor','k');
xlabel('Stimulus Intensity');
ylabel('proportion correct');
drawnow

%Function to be fitted
PF = @PAL_Logistic;

%Guesses for free parameters, fixed values for fixed parameters
params = [0 1 .5 0];    %or e.g.: [0 1 .5 0; 0 1 .5 0];

%Optional arguments for PAL_PFML_FitMultiple, 
%PAL_PFML_BootstrapParametricMultiple, 
%PAL_PFML_BootstrapNonParametricMultiple, PAL_PFLR_ModelComparison, and
%PAL_PFML_GoodnessOfFitMultiple
options = PAL_minimize('options');   %PAL_minimize search options
options.TolFun = 1e-12;     %Increase desired precision on LL
options.TolX = 1e-12;       %Increase desired precision on parameters
options.MaxFunEvals = 5000; %Allow more function evals
options.MaxIter = 5000;     %Allow more iterations
options.Display = 'off';    %suppress fminsearch messages
lapseLimits = [0 1];        %Range on lapse rates. Will go ignored here
                            %since lapse rate is not a free parameter
maxTries = 4;               %Try each fit at most four times        
rangeTries = [2 1.9 0 0];   %Range of random jitter to apply to initial 
                            %parameter values on retries of failed fits.

                            

%Fit lesser '1PF' model (constrained thresholds, constrained slopes, 
%   fixed guess rates, fixed lapse rates).
[paramsL LL exitflag output trash numParams1T1S] = PAL_PFML_FitMultiple(StimLevels, NumPos, ...
    OutOfNum, params, PF, 'thresholds','constrained','slopes',...
    'constrained','guessrates','fixed','lapserates','fixed',...
    'lapseLimits',lapseLimits,'SearchOptions',options);

%Fit fuller '2PF' model (unconstrained thresholds, unconstrained slopes, 
%   fixed guess rates, fixed lapse rates).
[paramsF LL exitflag output trash numParams2T2S] = PAL_PFML_FitMultiple(StimLevels, NumPos, ...
    OutOfNum, paramsL, PF, 'thresholds','unconstrained','slopes',...
    'unconstrained','guessrates','fixed','lapserates','fixed',...
    'lapseLimits',lapseLimits,'SearchOptions',options);

%plot fitted functions
ProportionCorrectModel = PF(paramsF(1,:),StimLevelsFineGrain);
plot(StimLevelsFineGrain,ProportionCorrectModel,'-','linewidth',2,...
    'color',[0 .7 0]);
ProportionCorrectModel = PF(paramsF(2,:),StimLevelsFineGrain);
plot(StimLevelsFineGrain,ProportionCorrectModel,'-','linewidth',2,...
    'color',[0 .7 0]);
xlabel('Stimulus Intensity');
ylabel('Proportion Correct');
set(h1, 'HandleVisibility','off');

%Threshold plot (SE bars to be added later)
plot(1:2,paramsF(1:2,1),'ko','markersize',10,'markerfacecolor','k');
h2 = gca;
axis([.5 2.5 -1 1]);
set(h2, 'units','pixels','position',[75 60 150 175]);
hold on;
set(h2, 'fontsize',16);
set(h2, 'Xtick',[1 2]);
xlabel('Condition');
ylabel('Threshold');
set(h2, 'HandleVisibility','off');

%Slope plot (SE bars to be added later)
plot(1:2,paramsF(1:2,2),'ko','markersize',10,'markerfacecolor','k');
h3 = gca;
axis([.5 2.5 0 2]);
set(h3, 'units','pixels','position',[300 60 150 175]);
hold on;
set(h3, 'fontsize',16);
set(h3, 'Xtick',[1 2]);
xlabel('Condition');
ylabel('Slope');
set(h3, 'HandleVisibility','off');
drawnow

%Determine Standard errors on Thresholds and Slopes
message = 'Determining standard errors......';
disp(message);
if ParOrNonPar == 1
    [SD paramsSim LLSim converged] = ...
        PAL_PFML_BootstrapParametricMultiple(StimLevels, OutOfNum, ...
        paramsF, Bse, PF, 'thresholds','unconstrained','slopes',...
        'unconstrained','guessrates','fixed','lapserates','fixed', ...
        'lapseLimits',lapseLimits,'SearchOptions',options,'maxTries',...
        maxTries,'rangeTries',rangeTries); 
else
    [SD paramsSim LLSim converged] = ...
        PAL_PFML_BootstrapNonParametricMultiple(StimLevels, NumPos, ...
        OutOfNum, paramsF, Bse, PF, 'thresholds','unconstrained',...
        'slopes','unconstrained','guessrates','fixed','lapserates',...
        'fixed', 'lapseLimits',lapseLimits,'SearchOptions',options,...
        'maxTries',maxTries,'rangeTries',rangeTries); 
end

%Add standard error bars to graphs
set(h2, 'HandleVisibility','on');
axes(h2);
line([1 1],[paramsF(1,1)-SD(1,1) paramsF(1,1)+SD(1,1)],'color','k',...
    'linewidth',2);
line([2 2],[paramsF(2,1)-SD(2,1) paramsF(2,1)+SD(2,1)],'color','k',...
    'linewidth',2);
set(h2, 'HandleVisibility','off');

set(h3, 'HandleVisibility','on');
axes(h3);
line([1 1],[paramsF(1,2)-SD(1,2) paramsF(1,2)+SD(1,2)],'color','k',...
    'linewidth',2);
line([2 2],[paramsF(2,2)-SD(2,2) paramsF(2,2)+SD(2,2)],'color','k',...
    'linewidth',2);
set(h2, 'HandleVisibility','off');
drawnow

message = 'Performing model comparison: Effect on threshold and/or slope?';
disp(message);

%omnibus test (2 PF vs 1 PF model, same threshold AND same slope?)
[TLR pTLR paramsL paramsF TLRSim converged] = ...
    PAL_PFLR_ModelComparison(StimLevels, NumPos, OutOfNum, params, Bmc, ...
    PF,'maxTries',maxTries, 'rangeTries',rangeTries,'lapseLimits', ...
    lapseLimits,'searchOptions',options);

%Plot fits under Fuller and Lesser models
figure('name','2 thresholds and 2 slopes vs. 1 threshold and 1 slope',...
    'units','pixels','position',[100 100 500 500]);
plot(StimLevels(1,:),ProportionCorrectObserved(1,:),'ko','markersize',...
    10,'markerfacecolor','k');
h1 = gca;
set(h1, 'units','pixels','position',[75 300 375 175]);
set(h1, 'fontsize',16);
set(h1, 'Xtick',StimLevels(1,:));
set(h1, 'Ytick',[.5:.1:1]);
axis([-2.5 2.5 .5 1]);
hold on;
plot(StimLevels(2,:),ProportionCorrectObserved(2,:),'ks','markersize',...
    10,'markerfacecolor','k');
ProportionCorrectModel = PF(paramsF(1,:),StimLevelsFineGrain);
plot(StimLevelsFineGrain,ProportionCorrectModel,'-','linewidth',2,...
    'color',[0 .7 0]);
ProportionCorrectModel = PF(paramsF(2,:),StimLevelsFineGrain);
plot(StimLevelsFineGrain,ProportionCorrectModel,'-','linewidth',2,...
    'color',[0 .7 0]);
ProportionCorrectModel = PF(paramsL(1,:),StimLevelsFineGrain);
plot(StimLevelsFineGrain,ProportionCorrectModel,'-','linewidth',2,...
    'color',[.7 0 0]);
xlabel('Stimulus Intensity');
ylabel('Proportion Correct');
text(-2.3,.95,'fuller','color',[0 .7 0],'FontSize',16);
text(-2.3,.9,'lesser','color',[.7 0 0],'FontSize',16);
set(gca,'HandleVisibility','off')

%Plot histogram of simulated TLR values and TLR value based on data. If
%function chi2pdf is detected add theoretical chi2 distribution with
%appropriate df.

[n centers] = hist(TLRSim,40);
hist(TLRSim,40)
h = findobj(gca,'Type','patch');
set(gca,'FOntSize',12)
set(h,'FaceColor','y','EdgeColor','k')
set(gca, 'units','pixels','position',[75 60 375 175]);
set(gca,'xlim',[0 1.2*max(TLR,centers(length(centers)))]);
xlim = get(gca, 'Xlim');
hold on
if exist('chi2pdf.m') == 2
    chi2x = xlim(1):xlim(2)/250:xlim(2);
    [maxim I]= max(n);
    chi2 = chi2pdf(chi2x,numParams2T2S-numParams1T1S)*(maxim/chi2pdf(centers(I),numParams2T2S-numParams1T1S));
    plot(chi2x,chi2,'k-','linewidth',2)
end
ylim = get(gca, 'Ylim');
plot(TLR,.05*ylim(2),'kv','MarkerSize',12,'MarkerFaceColor','k')
text(TLR,.15*ylim(2),'TLR data','Fontsize',11,'horizontalalignment',...
    'center');
message = ['p_{simul}: ' num2str(pTLR,'%5.4f')];
text(.95*xlim(2),.8*ylim(2),message,'horizontalalignment','right',...
    'fontsize',10);
if exist('chi2cdf.m') == 2
    message = ['p_{chi2}: ' num2str(1-chi2cdf(TLR,numParams2T2S-numParams1T1S),'%5.4f')];
    text(.95*xlim(2),.7*ylim(2),message,'horizontalalignment','right',...
        'fontsize',10);
end
xlabel('Simulated TLRs','FontSize',16)
ylabel('frequency','FontSize',16);
drawnow

message = 'Performing model comparison: Effect on threshold?';
disp(message);

%Fit lesser model (constrained thresholds, unconstrained slopes, 
%   fixed guess rates, fixed lapse rates).
[paramsL LL exitflag output trash numParams1T2S] = PAL_PFML_FitMultiple(StimLevels, NumPos, ...
    OutOfNum, params, PF, 'thresholds','constrained','slopes',...
    'unconstrained','guessrates','fixed','lapserates','fixed',...
    'lapseLimits',lapseLimits,'SearchOptions',options);

%thresholds (2 Thresholds, 2 Slopes vs 1 Threshold, 2 Slopes)
[TLR pTLR paramsL paramsF TLRSim converged] = ...
    PAL_PFLR_ModelComparison(StimLevels, NumPos, OutOfNum, params, Bmc, ...
    PF, 'lesserSlopes','unconstrained','maxTries',maxTries, ...
    'rangeTries',rangeTries,'lapseLimits', lapseLimits,'searchOptions',...
    options);

%Plot fits under Fuller and Lesser models
figure('name','2 thresholds and 2 slopes vs. 1 threshold and 2 slopes',...
    'units','pixels','position',[100 100 500 500]);
plot(StimLevels(1,:),ProportionCorrectObserved(1,:),'ko','markersize',...
    10,'markerfacecolor','k');
h1 = gca;
set(h1, 'units','pixels','position',[75 300 375 175]);
set(h1, 'fontsize',16);
set(h1, 'Xtick',StimLevels(1,:));
set(h1, 'Ytick',[.5:.1:1]);
axis([-2.5 2.5 .5 1]);
hold on;
plot(StimLevels(2,:),ProportionCorrectObserved(2,:),'ks','markersize',...
    10,'markerfacecolor','k');
ProportionCorrectModel = PF(paramsF(1,:),StimLevelsFineGrain);
plot(StimLevelsFineGrain,ProportionCorrectModel,'-','linewidth',2,...
    'color',[0 .7 0]);
ProportionCorrectModel = PF(paramsF(2,:),StimLevelsFineGrain);
plot(StimLevelsFineGrain,ProportionCorrectModel,'-','linewidth',2,...
    'color',[0 .7 0]);
ProportionCorrectModel = PF(paramsL(1,:),StimLevelsFineGrain);
plot(StimLevelsFineGrain,ProportionCorrectModel,'-','linewidth',2,...
    'color',[.7 0 0]);
ProportionCorrectModel = PF(paramsL(2,:),StimLevelsFineGrain);
plot(StimLevelsFineGrain,ProportionCorrectModel,'-','linewidth',2,...
    'color',[.7 0 0]);
xlabel('Stimulus Intensity');
ylabel('Proportion Correct');
text(-2.3,.95,'fuller','color',[0 .7 0],'FontSize',16);
text(-2.3,.9,'lesser','color',[.7 0 0],'FontSize',16);
set(gca,'HandleVisibility','off')

%Plot histogram of simulated TLR values, TLR value based on data. If
%function chi2pdf is detected add theoretical chi2 distribution with
%appropriate df.

[n centers] = hist(TLRSim,40);
hist(TLRSim,40)
h = findobj(gca,'Type','patch');
set(gca,'FOntSize',12)
set(h,'FaceColor','y','EdgeColor','k')
set(gca, 'units','pixels','position',[75 60 375 175]);
set(gca,'xlim',[0 1.2*max(TLR,centers(length(centers)))]);
xlim = get(gca, 'Xlim');
hold on
if exist('chi2pdf.m') == 2
    chi2x = xlim(1):xlim(2)/250:xlim(2);
    [maxim I]= max(n);
    chi2 = chi2pdf(chi2x,numParams2T2S-numParams1T2S)*(n(2)/chi2pdf(centers(2),numParams2T2S-numParams1T2S));
    plot(chi2x,chi2,'k-','linewidth',2)
end
ylim = get(gca, 'Ylim');
plot(TLR,.05*ylim(2),'kv','MarkerSize',12,'MarkerFaceColor','k')
text(TLR,.15*ylim(2),'TLR data','Fontsize',11,'horizontalalignment',...
    'center');
message = ['p_{simul}: ' num2str(pTLR,'%5.4f')];
text(.95*xlim(2),.8*ylim(2),message,'horizontalalignment','right',...
    'fontsize',10);
if exist('chi2cdf.m') == 2
    message = ['p_{chi2}: ' num2str(1-chi2cdf(TLR,numParams2T2S-numParams1T2S),'%5.4f')];
    text(.95*xlim(2),.7*ylim(2),message,'horizontalalignment','right',...
        'fontsize',10);
end    
xlabel('Simulated TLRs','FontSize',16)
ylabel('frequency','FontSize',16);
drawnow

message = 'Performing model comparison: Effect on slope?';
disp(message);

%Fit lesser model (unconstrained thresholds, constrained slopes, 
%   fixed guess rates, fixed lapse rates).
[paramsL LL exitflag output trash numParams2T1S] = PAL_PFML_FitMultiple(StimLevels, NumPos, ...
    OutOfNum, params, PF, 'thresholds','unconstrained','slopes',...
    'constrained','guessrates','fixed','lapserates','fixed',...
    'lapseLimits',lapseLimits,'SearchOptions',options);

%slopes (2 Thresholds, 2 Slopes vs 2 Thresholds, 1 Slope)
[TLR pTLR paramsL paramsF TLRSim converged] = ...
    PAL_PFLR_ModelComparison(StimLevels, NumPos, OutOfNum, params, Bmc, ...
    PF, 'lesserThresholds','unconstrained','maxTries',maxTries, ...
    'rangeTries',rangeTries,'lapseLimits', lapseLimits,'searchOptions',...
    options);

%Plot fits under Fuller and Lesser models
figure('name','2 thresholds and 2 slopes vs. 2 threshold and 1 slopes',...
    'units','pixels','position',[100 100 500 500]);
plot(StimLevels(1,:),ProportionCorrectObserved(1,:),'ko','markersize',...
    10,'markerfacecolor','k');
h1 = gca;
set(h1, 'units','pixels','position',[75 300 375 175]);
set(h1, 'fontsize',16);
set(h1, 'Xtick',StimLevels(1,:));
set(h1, 'Ytick',[.5:.1:1]);
axis([-2.5 2.5 .5 1]);
hold on;
plot(StimLevels(2,:),ProportionCorrectObserved(2,:),'ks','markersize',...
    10,'markerfacecolor','k');
ProportionCorrectModel = PF(paramsF(1,:),StimLevelsFineGrain);
plot(StimLevelsFineGrain,ProportionCorrectModel,'-','linewidth',2,...
    'color',[0 .7 0]);
ProportionCorrectModel = PF(paramsF(2,:),StimLevelsFineGrain);
plot(StimLevelsFineGrain,ProportionCorrectModel,'-','linewidth',2,...
    'color',[0 .7 0]);
ProportionCorrectModel = PF(paramsL(1,:),StimLevelsFineGrain);
plot(StimLevelsFineGrain,ProportionCorrectModel,'--','linewidth',2,...
    'color',[.7 0 0]);
ProportionCorrectModel = PF(paramsL(2,:),StimLevelsFineGrain);
plot(StimLevelsFineGrain,ProportionCorrectModel,'--','linewidth',2,...
    'color',[.7 0 0]);
xlabel('Stimulus Intensity');
ylabel('Proportion Correct');
text(-2.3,.95,'fuller','color',[0 .7 0],'FontSize',16);
text(-2.3,.9,'lesser','color',[.7 0 0],'FontSize',16);
set(gca,'HandleVisibility','off')

%Plot histogram of simulated TLR values, TLR value based on data. If
%function chi2pdf is detected add theoretical chi2 distribution with
%appropriate df.

[n centers] = hist(TLRSim,40);
hist(TLRSim,40)
h = findobj(gca,'Type','patch');
set(gca,'FOntSize',12)
set(h,'FaceColor','y','EdgeColor','k')
set(gca, 'units','pixels','position',[75 60 375 175]);
set(gca,'xlim',[0 1.2*max(TLR,centers(length(centers)))]);
xlim = get(gca, 'Xlim');
hold on
if exist('chi2pdf.m') == 2
    chi2x = xlim(1):xlim(2)/250:xlim(2);
    [maxim I]= max(n);
    chi2 = chi2pdf(chi2x,numParams2T2S-numParams2T1S)*(n(2)/chi2pdf(centers(2),numParams2T2S-numParams2T1S));
    plot(chi2x,chi2,'k-','linewidth',2)
end
ylim = get(gca, 'Ylim');
plot(TLR,.05*ylim(2),'kv','MarkerSize',12,'MarkerFaceColor','k')
text(TLR,.15*ylim(2),'TLR data','Fontsize',11,'horizontalalignment',...
    'center');
message = ['p_{simul}: ' num2str(pTLR,'%5.4f')];
text(.95*xlim(2),.8*ylim(2),message,'horizontalalignment','right',...
    'fontsize',10);
if exist('chi2cdf.m') == 2
    message = ['p_{chi2}: ' num2str(1-chi2cdf(TLR,numParams2T2S-numParams2T1S),'%5.4f')];
    text(.95*xlim(2),.7*ylim(2),message,'horizontalalignment','right',...
        'fontsize',10);
end
xlabel('Simulated TLRs','FontSize',16)
ylabel('frequency','FontSize',16);
drawnow


message = sprintf('Performing model comparison: 2 thresholds, 1 slope v');
message = strcat(message,'s. saturated model (i.e., Goodness-of-Fit)');
disp(message);

numParamsSat = sum(sum(OutOfNum~=0));

%Goodness-of-fit
[TLR pTLR TLRSim converged] =PAL_PFML_GoodnessOfFitMultiple(StimLevels, ...
    NumPos, OutOfNum, paramsL, Bmc, PF, 'Thresholds', 'unconstrained', ...
    'Slopes', 'constrained', 'GuessRates', 'fixed', 'LapseRates', ...
    'fixed','maxTries',maxTries, 'rangeTries',rangeTries,'lapseLimits', ...
    lapseLimits,'searchOptions', options);


%Plot fits under Fuller and Lesser models
figure('name','Saturated vs. 2 threshold and 1 slopes','units','pixels',...
    'position',[100 100 500 500]);
plot(StimLevels(1,:),ProportionCorrectObserved(1,:),'o','color',...
    [0 .7 0],'markersize',10,'markerfacecolor',[0 .7 0]);
h1 = gca;
set(h1, 'units','pixels','position',[75 300 375 175]);
set(h1, 'fontsize',16);
set(h1, 'Xtick',StimLevels(1,:));
set(h1, 'Ytick',[.5:.1:1]);
axis([-2.5 2.5 .5 1]);
hold on;
plot(StimLevels(2,:),ProportionCorrectObserved(2,:),'s','color',...
    [0 .7 0],'markersize',10,'markerfacecolor',[0 .7 0]);
ProportionCorrectModel = PF(paramsL(1,:),StimLevelsFineGrain);
plot(StimLevelsFineGrain,ProportionCorrectModel,'-','linewidth',2,...
    'color',[.7 0 0]);
ProportionCorrectModel = PF(paramsL(2,:),StimLevelsFineGrain);
plot(StimLevelsFineGrain,ProportionCorrectModel,'-','linewidth',2,...
    'color',[.7 0 0]);
xlabel('Stimulus Intensity');
ylabel('Proportion Correct');
text(-2.3,.95,'fuller','color',[0 .7 0],'FontSize',16);
text(-2.3,.9,'lesser','color',[.7 0 0],'FontSize',16);
set(gca,'HandleVisibility','off')

%Plot histogram of simulated TLR values, TLR value based on data. If
%function chi2pdf is detected add theoretical chi2 distribution with
%appropriate df.
[n centers] = hist(TLRSim,40);
hist(TLRSim,40)
set(gca,'FOntSize',12)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','y','EdgeColor','k')
set(gca, 'units','pixels','position',[75 60 375 175]);
set(gca,'xlim',[0 1.2*max(TLR,centers(length(centers)))]);
xlim = get(gca, 'Xlim');
hold on
if exist('chi2pdf.m') == 2
    chi2x = xlim(1):xlim(2)/250:xlim(2);
    [maxim I]= max(n);
    chi2 = chi2pdf(chi2x,numParamsSat-numParams2T1S)*(maxim/chi2pdf(centers(I),numParamsSat-numParams2T1S));
    plot(chi2x,chi2,'k-','linewidth',2)
end
ylim = get(gca, 'Ylim');
plot(TLR,.05*ylim(2),'kv','MarkerSize',12,'MarkerFaceColor','k')
text(TLR,.15*ylim(2),'TLR data','Fontsize',11,'horizontalalignment',...
    'center');
message = ['p_{simul}: ' num2str(pTLR,'%5.4f')];
text(.95*xlim(2),.8*ylim(2),message,'horizontalalignment','right',...
    'fontsize',10);
if exist('chi2cdf.m') == 2
    message = ['p_{chi2}: ' num2str(1-chi2cdf(TLR,numParamsSat-numParams2T1S),'%5.4f')];
    text(.95*xlim(2),.7*ylim(2),message,'horizontalalignment','right',...
        'fontsize',10);
end
xlabel('Simulated TLRs','FontSize',16)
ylabel('frequency','FontSize',16);

drawnow
toc