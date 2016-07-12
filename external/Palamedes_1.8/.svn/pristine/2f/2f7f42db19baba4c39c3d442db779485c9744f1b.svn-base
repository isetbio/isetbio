%
%PAL_PFLR_FourGroup_Demo  Demonstrates use of Palamedes functions to 
%(1) fit psychometric functions to multiple conditions simultaneously
%while allowing considerable flexibility to constrain parameter values
%between conditions (including the use of contrasts to do so)(2) determine
%standard errors of the free parameters in simultaneous fits using
%bootstrap simulations and (3) perform statistical model comparisons 
%between nested models. Different model comparisons are performed including
%a goodness-of-fit test for a user-defined model.
%
%Program first fits psychometric functions (PFs) to some data in a four
%group experiment. The model fits individual thresholds to the four
%conditions, but assumes that the slopes and lapse rates of the PFs are 
%identical between conditions. It assumes the guess rate is 0.5. Thus, this
%is a 6 parameter model (4 thresholds, 1 slope, 1 lapse rate). Plots of the 
%fit are presented in a figure.
%
%Program then estimates standard error on the 6 free parameters using
%bootstrap simulations. Standard error bars are added to the figure.
%
%Next a model comparison is performed. The model above is compared to a
%model which assumes thresholds are identical between conditions. A figure
%is created which displays parameter estimates for both models. Also
%displayed is the distribution of simulated transformed likelihood ratios.
%[TLR = -2x(log likelihood lesser model - log likelihood fuller model)].
%
%A second model comparison tests whether the decelerating trend in the
%threshold plot is real or whether a linear trend suffices. A figure is
%created similar to that for the previous model comparison.
%
%Finally, a model which assumes that thresholds follow a quadratic trend,
%that slopes and lapse rates are equal between conditions is compared
%against the saturated model (i.e., this is a Goodness-of-fit test). Figure
%similar to those above is created.
%
%For all model comparisons, the proportion of simulated TLR values greater 
%than the data's TLR value is displayed.
%
%In case the Matlab routines chi2cdf and chi2pdf are present, program
%superimposes theoretical chi-square distribution with appropriate degrees 
%of freedom on empirical sampling distribution and presents the p-value 
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
%PAL_Contrasts
%PAL_Logistic
%
%More information on any of these functions may be found by typing
%help followed by the name of the function. e.g., help PAL_PFML_FitMultiple
%
%NP (September 2009)

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

%Stimulus Levels 4 (conditions) x 5 (stimulus levels)
%rows need not be identical
StimLevels = repmat([-2:1:2],[4 1]);

%Number of trials at each entry of StimLevels
%entries may be zero (useful if stimulus levels are not identical between
%conditions).
OutOfNum = 150*ones(4,5);

%Number of positive responses (e.g., 'yes', 'correct') for each entry of 
%StimLevels
NumPos =   [84 92 128 137 143; ...
            67 85 103 131 139; ...
            73 85  92 125 143; ...
            82 86  97 122 141];

%Plot raw data
ProportionCorrectObserved = NumPos ./ OutOfNum; 
StimLevelsFineGrain = [min(min(StimLevels)):(max(max(StimLevels) - ... 
    min(min(StimLevels))))./1000:max(max(StimLevels))];
figure('name','Model: 4 thresholds, shared slope and lapse rate',...
    'units','pixels','position',[100 100 500 500]);
plot(StimLevels(1,:),ProportionCorrectObserved(1,:),'bo','markersize',...
    10,'markerfacecolor','b');
hold on
plot(StimLevels(2,:),ProportionCorrectObserved(2,:),'ro','markersize',...
    10,'markerfacecolor','r');
plot(StimLevels(3,:),ProportionCorrectObserved(3,:),'go','markersize',...
    10,'markerfacecolor','g');
plot(StimLevels(4,:),ProportionCorrectObserved(4,:),'mo','markersize',...
    10,'markerfacecolor','m');
h1 = gca;
set(h1, 'units','pixels','position',[75 300 375 175]);
set(h1, 'fontsize',16);
set(h1, 'Xtick',StimLevels(1,:));
set(h1, 'Ytick',[.5:.1:1]);
axis([-2.5 2.5 .5 1]);
hold on;
xlabel('Stimulus Intensity');
ylabel('proportion correct');
drawnow

%Function to be fitted
PF = @PAL_Logistic;

%Guesses for free parameters, fixed values for fixed parameters
params = [0 2 .5 0.02];       

%or, e.g.,
%   params = [-.6 1.8 .5 .02; ...
%              .1 1.8 .5 .02; ...
%              .6 1.8 .5 .02; ...
%              .9 1.8 .5 .02];


           
%Optional arguments for PAL_PFML_FitMultiple, 
%PAL_PFML_BootstrapParametricMultiple, 
%PAL_PFML_BootstrapNonParametricMultiple, PAL_PFLR_ModelComparison, and
%PAL_PFML_GoodnessOfFitMultiple
options = PAL_minimize('options');   %Nelder-Mead search options
options.MaxFunEvals = 5000;         %Allow more function evals
options.MaxIter = 5000;             %Allow more iterations
options.TolFun = 1e-12;             %Increase desired precision on LL
options.TolX = 1e-12;               %Increase desired precision on params
                                    %default precision (1e-04) results in
                                    %poor estimates here.
options.Display = 'notify';            %suppress fminsearch messages

lapseLimits = [0 .5];         %Range on lapse rates.  
maxTries = 10;                %Try each fit at most ten times        
rangeTries = [2 2 0 0.04];    %Range of random jitter to apply to initial 
                              %parameter values on retries of failed fits.


%Fit model that assumes slopes are equal between conditions, guessrates are
%0.5, lapse rates are equal between conditions, but allows thresholds to
%differ between conditions. 6 free parameters: 4 thresholds, 1 slope, 1
%lapse rate.
[params4T1S1L LL exitflag output trash numParams4T1S1L] = PAL_PFML_FitMultiple(StimLevels, ...
    NumPos, OutOfNum, params, PF, 'Thresholds', 'unconstrained', ...
    'Slopes', 'constrained', 'GuessRates', 'fixed', 'LapseRates', ...
    'constrained','lapseLimits', lapseLimits,'SearchOptions',options);

%Add fitted psychometric functions to raw data plot
ProportionCorrectModel = PF(params4T1S1L(1,:),StimLevelsFineGrain);
plot(StimLevelsFineGrain,ProportionCorrectModel,'-','linewidth',2,...
    'color','b');
ProportionCorrectModel = PF(params4T1S1L(2,:),StimLevelsFineGrain);
plot(StimLevelsFineGrain,ProportionCorrectModel,'-','linewidth',2,...
    'color','r');
ProportionCorrectModel = PF(params4T1S1L(3,:),StimLevelsFineGrain);
plot(StimLevelsFineGrain,ProportionCorrectModel,'-','linewidth',2,...
    'color','g');
ProportionCorrectModel = PF(params4T1S1L(4,:),StimLevelsFineGrain);
plot(StimLevelsFineGrain,ProportionCorrectModel,'-','linewidth',2,...
    'color','m');
set(h1, 'HandleVisibility','off');

%Threshold plot (SE bars to be added later)
plot(1,params4T1S1L(1,1),'bo','markersize',10,'markerfacecolor','b');
hold on;
plot(2,params4T1S1L(2,1),'ro','markersize',10,'markerfacecolor','r');
plot(3,params4T1S1L(3,1),'go','markersize',10,'markerfacecolor','g');
plot(4,params4T1S1L(4,1),'mo','markersize',10,'markerfacecolor','m');
h2 = gca;
axis([.5 4.5 -1 1]);
set(h2, 'units','pixels','position',[75 60 175 175]);
hold on;
set(h2, 'fontsize',12);
set(h2, 'Xtick',[1 2 3 4]);
xlabel('Condition');
ylabel('Threshold');
set(h2, 'HandleVisibility','off');
drawnow

%Slope plot (SE bar to be added later)
plot(1,params4T1S1L(1,2),'ko','markersize',10,'markerfacecolor','k');
h3 = gca;
axis([.5 1.5 1 3]);
set(h3, 'units','pixels','position',[300 60 50 175]);
hold on;
set(h3, 'fontsize',12);
set(h3, 'Xtick',[],'ytick',[1 2 3]);
xlabel('All Conditions');
ylabel('Slope');
set(h3, 'HandleVisibility','off');
drawnow

%Lapse Rate plot (SE bar to be added later)
plot(1,params4T1S1L(1,4),'ko','markersize',10,'markerfacecolor','k');
h4 = gca;
axis([.5 1.5 0 .08]);
set(h4, 'units','pixels','position',[400 60 50 175]);
hold on;
set(h4, 'fontsize',12);
set(h4, 'Xtick',[],'Ytick',[0:.02:.08],'YtickLabel',{'0','2','4','6','8'});
text(.5,.088,'x 10^{-2}');
xlabel('All Conditions');
ylabel('Lapse Rate');
set(h4, 'HandleVisibility','off');
drawnow

%Determine Standard errors on free parameters

message = 'Determining standard errors......';
disp(message);
tic
if ParOrNonPar == 1
    [SD paramsSim LLSim converged] = ...
        PAL_PFML_BootstrapParametricMultiple(StimLevels, OutOfNum, ...
        params4T1S1L, Bse, PF, 'Thresholds', 'unconstrained', 'Slopes', ...
        'constrained', 'GuessRates', 'fixed', 'LapseRates', ...
        'constrained','maxTries',maxTries, 'rangeTries',rangeTries,...
        'lapseLimits', lapseLimits,'SearchOptions',options);
else
    [SD paramsSim LLSim converged] = ...
        PAL_PFML_BootstrapNonParametricMultiple(StimLevels, NumPos, ...
        OutOfNum, params4T1S1L, Bse, PF, 'Thresholds', 'unconstrained', ...
        'Slopes', 'constrained', 'GuessRates', 'fixed', 'LapseRates', ...
        'constrained','maxTries',maxTries, 'rangeTries',rangeTries,...
        'lapseLimits', lapseLimits,'SearchOptions',options);
end
toc

%Add standard error bars to graphs
%Threshold standard errors
set(h2, 'HandleVisibility','on');
axes(h2);
line([1 1],[params4T1S1L(1,1)-SD(1,1) params4T1S1L(1,1)+SD(1,1)],...
    'color','b','linewidth',2);
line([2 2],[params4T1S1L(2,1)-SD(2,1) params4T1S1L(2,1)+SD(2,1)],...
    'color','r','linewidth',2);
line([3 3],[params4T1S1L(3,1)-SD(3,1) params4T1S1L(3,1)+SD(3,1)],...
    'color','g','linewidth',2);
line([4 4],[params4T1S1L(4,1)-SD(4,1) params4T1S1L(4,1)+SD(4,1)],...
    'color','m','linewidth',2);
set(h2, 'HandleVisibility','off');

%Slope standard error
set(h3, 'HandleVisibility','on');
axes(h3);
line([1 1],[params4T1S1L(1,2)-SD(1,2) params4T1S1L(1,2)+SD(1,2)],...
    'color','k','linewidth',2);
set(h2, 'HandleVisibility','off');

%Lapse rate standard errors
set(h4, 'HandleVisibility','on');
axes(h4);
line([1 1],[params4T1S1L(1,4)-SD(1,4) params4T1S1L(1,4)+SD(1,4)],...
    'color','k','linewidth',2);
set(h2, 'HandleVisibility','off');
drawnow

%Test whether thresholds differ (i.e., test 6 parameter model (4
%thresholds, 1 slope, 1 lapse rate) against 3 parameter model (1 threshold, 
%1 slope, 1 lapse rate) model.

message = 'Performing model comparison: thresholds identical?.....';
disp(message);

[params1T1S1L LL exitflag output trash numParams1T1S1L] = PAL_PFML_FitMultiple(StimLevels, ...
    NumPos, OutOfNum, params, PF, 'Thresholds', 'constrained', ...
    'Slopes', 'constrained', 'GuessRates', 'fixed', 'LapseRates', ...
    'constrained','lapseLimits', lapseLimits,'SearchOptions',options);

[TLR pTLR params1T1S1L params4T1S1L TLRSim converged] = ...
    PAL_PFLR_ModelComparison (StimLevels, NumPos, OutOfNum, ...
    params4T1S1L, Bmc, PF, 'lesserlapse','constrained', 'fullerlapse',...
    'constrained', 'fullerslope','constrained','maxTries',maxTries, ...
    'rangeTries',rangeTries,'searchoptions',options,'lapseLimits', ...
    lapseLimits);

%Create figure that shows parameter estimates for both models, distribution 
%of simulated TLR (transformed likelihood ratio), p-value based on 
%simulations and, if chi2pdf and chi2cdf are detected, the theoretical 
%chi-square distribution with appropriate df and p-value based on 
%theoretical chi-square.
figure('name','Thresholds differ?','units','pixels','position',...
    [100 100 500 500]);

%Threshold plot
plot(1:4,params4T1S1L(1:4,1),'o','color',[0 .7 0],'markersize',14,...
    'linewidth',2);
hold on;
plot(1:4,params1T1S1L(1:4,1),'o','color',[.7 0 0],'markersize',14,...
    'linewidth',2);
plot(1:4,params4T1S1L(1:4,1),'ko','markersize',10,'markerFaceColor','k');
h2 = gca;
axis([.5 4.5 -1 1]);
set(h2, 'units','pixels','position',[75 300 175 175]);
hold on;
set(h2, 'fontsize',12);
set(h2, 'Xtick',[1 2 3 4]);
xlabel('Condition');
ylabel('Threshold');
text(.7,.8,'fuller','color',[0 .7 0],'fontsize',12);
text(.7,.6,'lesser','color',[.7 0 0],'fontsize',12);
text(2,-.8,'unconstrained','color',[0 0 0],'fontsize',12);
set(h2, 'HandleVisibility','off');

%Slope plot
plot(1,params4T1S1L(1,2),'o','color',[0 .7 0],'markersize',14,'linewidth',2);
hold on
plot(1,params1T1S1L(1,2),'o','color',[.7 0 0],'markersize',14,'linewidth',2);
plot(1,params4T1S1L(1,2),'ko','markersize',10,'markerfacecolor','k');
h3 = gca;
axis([.5 1.5 1 3]);
set(h3, 'units','pixels','position',[300 300 50 175]);
hold on;
set(h3, 'fontsize',12);
set(h3, 'Xtick',[],'ytick',[1 2 3]);
xlabel('All Conditions');
ylabel('Slope');
set(h3, 'HandleVisibility','off');

%Lapse Rate plot
plot(1,params4T1S1L(1,4),'o','color',[0 .7 0],'markersize',14,'linewidth',2);
hold on
plot(1,params1T1S1L(1,4),'o','color',[.7 0 0],'markersize',14,'linewidth',2);
plot(1,params4T1S1L(1,4),'ko','markersize',10,'markerfacecolor','k');
h4 = gca;
axis([.5 1.5 0 .08]);
set(h4, 'units','pixels','position',[400 300 50 175]);
hold on;
set(h4, 'fontsize',12);
set(h4, 'Xtick',[],'Ytick',[0:.02:.08],'YtickLabel',{'0','2','4','6','8'});
text(.5,.088,'x 10^{-2}');
xlabel('All Conditions');
ylabel('Lapse Rate');
set(h4, 'HandleVisibility','off');

%Distribution of simulated TLRs
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
    chi2 = chi2pdf(chi2x,numParams4T1S1L-numParams1T1S1L)*...
        (maxim/chi2pdf(centers(I),numParams4T1S1L-numParams1T1S1L));
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
    message = ['p_{chi2}: ' num2str(1-chi2cdf(TLR,numParams4T1S1L-...
        numParams1T1S1L),'%5.4f')];
    text(.95*xlim(2),.7*ylim(2),message,'horizontalalignment','right',...
        'fontsize',10);
end
xlabel('Simulated TLRs','FontSize',16)
ylabel('frequency','FontSize',16);
drawnow

%Compare model assuming quadratic trend of thresholds against model
%assuming linear trend. I.e., is the bend apparent in threshold plot real?

message = 'Performing model comparison: Quadratic trend on thresholds?...';
disp(message);

Contrasts = PAL_Contrasts(4, 'polynomial'); %generate polynomial contrasts
ContrastsModelC = Contrasts(1:2,:);         %intercept + linear trend

%Fit model assuming linear trend (not necessary: PAL_PFLR_ModelComparison
%will perform fit. Included for demonstration only and to determine number
%of free parameters in model
[paramsC LLC exitflag output trash numParamsC] = ...
     PAL_PFML_FitMultiple(StimLevels, NumPos, OutOfNum, params, PF,...
    'Thresholds',ContrastsModelC,'Slopes','constrained','GuessRates',...
    'fixed', 'LapseRates', 'constrained','lapseLimits', lapseLimits,...
    'SearchOptions',options);

ContrastsModelD = Contrasts(1:3,:);         %intercept + linear + quadratic

%Fit model assuming quadratic trend (not necessary: PAL_PFLR_ModelComparison
%will perform fit. Included for demonstration only and to determine number
%of free parameters in model
[paramsD LLD exitflag output trash numParamsD] = ...
    PAL_PFML_FitMultiple(StimLevels, NumPos, OutOfNum, params, PF,...
    'Thresholds',ContrastsModelD, 'Slopes', 'constrained', 'GuessRates',...
    'fixed', 'LapseRates', 'constrained','lapseLimits', lapseLimits,...
    'SearchOptions',options);

%Perform actual model comparison

[TLR pTLR paramsC paramsD TLRSim converged] = ...
    PAL_PFLR_ModelComparison(StimLevels, NumPos, OutOfNum, paramsD, Bmc,...
    PF, 'lesserthreshold', ContrastsModelC, 'fullerthreshold', ...
    ContrastsModelD, 'lesserlapse','constrained', 'fullerlapse',...
    'constrained', 'fullerslope','constrained','maxTries',maxTries,...
    'rangeTries',rangeTries,'lapseLimits', lapseLimits,'SearchOptions',...
    options);

%Create figure that shows parameter estimates for both models, distribution 
%of simulated TLR (transformed likelihood ratio), p-value based on 
%simulations and, if chi2pdf and chi2cdf are detected, the theoretical 
%chi-square distribution with appropriate df and p-value based on 
%theoretical chi-square.
figure('name','Quadratic vs. Linear trend on thresholds (is bend real?)',...
    'units','pixels','position',[100 100 500 500]);

%Threshold plot
plot(1:4,paramsD(1:4,1),'o','color',[0 .7 0],'markersize',14,'linewidth',2);
hold on;
plot(1:4,paramsC(1:4,1),'o','color',[.7 0 0],'markersize',14,'linewidth',2);
plot(1:4,params4T1S1L(1:4,1),'ko','markersize',10,'markerFacecolor','k');
h2 = gca;
axis([.5 4.5 -1 1]);
set(h2, 'units','pixels','position',[75 300 175 175]);
hold on;
set(h2, 'fontsize',12);
set(h2, 'Xtick',[1 2 3 4]);
xlabel('Condition');
ylabel('Threshold');
text(.7,.8,'fuller','color',[0 .7 0],'fontsize',12);
text(.7,.6,'lesser','color',[.7 0 0],'fontsize',12);
text(2,-.8,'unconstrained','color',[0 0 0],'fontsize',12);
set(h2, 'HandleVisibility','off');

%Slope plot
plot(1,paramsD(1,2),'o','color',[0 .7 0],'markersize',14,'linewidth',2);
hold on
plot(1,paramsC(1,2),'o','color',[.7 0 0],'markersize',14,'linewidth',2);
plot(1,params4T1S1L(1,2),'ko','markersize',10,'markerfacecolor','k');
h3 = gca;
axis([.5 1.5 1 3]);
set(h3, 'units','pixels','position',[300 300 50 175]);
hold on;
set(h3, 'fontsize',12);
set(h3, 'Xtick',[],'ytick',[1 2 3]);
xlabel('All Conditions');
ylabel('Slope');
set(h3, 'HandleVisibility','off');

%Lapse Rate plot
plot(1,paramsD(1,4),'o','color',[0 .7 0],'markersize',14,'linewidth',2);
hold on
plot(1,paramsC(1,4),'o','color',[.7 0 0],'markersize',14,'linewidth',2);
plot(1,params4T1S1L(1,4),'ko','markersize',10,'markerfacecolor','k');
h4 = gca;
axis([.5 1.5 0 .08]);
set(h4, 'units','pixels','position',[400 300 50 175]);
hold on;
set(h4, 'fontsize',12);
set(h4, 'Xtick',[],'Ytick',[0:.02:.08],'YtickLabel',{'0','2','4','6','8'});
text(.5,.088,'x 10^{-2}');
xlabel('All Conditions');
ylabel('Lapse Rate');
set(h4, 'HandleVisibility','off');

%Distribution of simulated TLRs
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
    chi2 = chi2pdf(chi2x,numParamsD-numParamsC)*(n(2)/chi2pdf(centers(2),...
        numParamsD-numParamsC));
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
    message = ['p_{chi2}: ' num2str(1-chi2cdf(TLR,numParamsD-numParamsC),...
        '%5.4f')];
    text(.95*xlim(2),.7*ylim(2),message,'horizontalalignment','right',...
        'fontsize',10);
end
xlabel('Simulated TLRs','FontSize',16)
ylabel('frequency','FontSize',16);
drawnow

%Goodness of fit of quadratic model.

message = 'Performing model comparison: Goodness-of-fit of quadratic model...';
disp(message);

[TLR pTLR TLRSim converged] = PAL_PFML_GoodnessOfFitMultiple(StimLevels,...
    NumPos, OutOfNum, paramsD, Bmc, PF, 'Thresholds', ContrastsModelD,...
    'Slopes', 'constrained', 'GuessRates', 'fixed', 'LapseRates',...
    'constrained', 'maxTries', maxTries, 'rangeTries',rangeTries,...
    'lapseLimits', lapseLimits,'SearchOptions',options);

numParamsSat = sum(sum(OutOfNum~=0));

%Create figure that shows parameter estimates for both models, distribution 
%of simulated TLR (transformed likelihood ratio), p-value based on 
%simulations and, if chi2pdf and chi2cdf are detected, the theoretical 
%chi-square distribution with appropriate df and p-value based on 
%theoretical chi-square.
figure('name','Goodness-of-fit','units','pixels',...
    'position',[100 100 500 500]);
plot(StimLevels(1,:),ProportionCorrectObserved(1,:),'o','color',...
    [0 .7 0],'markersize',10,'markerfacecolor',[0 .7 0]);
hold on
plot(StimLevels(2,:),ProportionCorrectObserved(2,:),'o','color',...
    [0 .7 0],'markersize',10,'markerfacecolor',[0 .7 0]);
plot(StimLevels(3,:),ProportionCorrectObserved(3,:),'o','color',...
    [0 .7 0],'markersize',10,'markerfacecolor',[0 .7 0]);
plot(StimLevels(4,:),ProportionCorrectObserved(4,:),'o','color',...
    [0 .7 0],'markersize',10,'markerfacecolor',[0 .7 0]);
h1 = gca;
set(h1, 'units','pixels','position',[75 300 375 175]);
set(h1, 'fontsize',16);
set(h1, 'Xtick',StimLevels(1,:));
set(h1, 'Ytick',[.5:.1:1]);
axis([-2.5 2.5 .5 1]);
hold on;
xlabel('Stimulus Intensity');
ylabel('proportion correct');
text(-2.25,.95,'fuller','color',[0 .7 0],'fontsize',12);
text(-2.25,.9,'lesser','color',[.7 0 0],'fontsize',12);

ProportionCorrectModel = PF(paramsD(1,:),StimLevelsFineGrain);
plot(StimLevelsFineGrain,ProportionCorrectModel,'-','linewidth',2,...
    'color',[.7 0 0]);
ProportionCorrectModel = PF(paramsD(2,:),StimLevelsFineGrain);
plot(StimLevelsFineGrain,ProportionCorrectModel,'-','linewidth',2,...
    'color',[.7 0 0]);
ProportionCorrectModel = PF(paramsD(3,:),StimLevelsFineGrain);
plot(StimLevelsFineGrain,ProportionCorrectModel,'-','linewidth',2,...
    'color',[.7 0 0]);
ProportionCorrectModel = PF(paramsD(4,:),StimLevelsFineGrain);
plot(StimLevelsFineGrain,ProportionCorrectModel,'-','linewidth',2,...
    'color',[.7 0 0]);
set(h1, 'HandleVisibility','off');

%Distribution of simulated TLRs
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
    chi2 = chi2pdf(chi2x,numParamsSat-numParamsD)*...
        (maxim/chi2pdf(centers(I),numParamsSat-numParamsD));
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
    message = ['p_{chi2}: ' num2str(1-chi2cdf(TLR,numParamsSat-...
        numParamsD),'%5.4f')];
    text(.95*xlim(2),.7*ylim(2),message,'horizontalalignment','right',...
        'fontsize',10);
end
xlabel('Simulated TLRs','FontSize',16)
ylabel('frequency','FontSize',16);
drawnow