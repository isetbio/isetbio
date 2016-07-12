%
%PAL_PFML_gammaEQlambda_Demo demonstrates the option to constrain the guess
%rate and the lapse rate to be equal. It does so in the context of the fit
%to a single PF. Usage generalizes to fitting multiple PFs simultaneously
%using PAL_PFML_FitMultiple (and associated functions:
%PAL_PFML_BootstrapParametricMultiple,
%PAL_PFML_BootstrapNonParametricMultiple,
%PAL_PFML_GoodnessOfFitMultiple
%PAL_PFLR_ModelComparison).
%
%The option to constrain gamma and lambda parameters to be equal was 
%introduced in Palamedes version 1.3.0.
%
%Demonstrates usage of Palamedes functions:
%-PAL_PFML_Fit
%-PAL_PFML_BootstrapParametric
%-PAL_PFML_BootstrapNonParametric
%-PAL_PFML_GoodnessOfFit
%-PAL_PFML_FitMultiple
%The above are demonstrated here especially with regard to usage of the
%optional function argument 'gammaEQlambda'.
%
%NP (September 2011)

clear all

rand('state',sum(100*clock));

disp([sprintf('\n')  'A Psychometric Function model will be fit in several ways: (1) All 4']);
disp(['parameters are free to vary [green curve] (2)guess and lapse rate are']);
disp(['constrained to be equal [red curve] (3) guess and lapse rate are constrained']);
disp(['to be equal and ''iAPLE'' option to incorporate lapse rates is used [blue']);
disp(['curve]. Curves are in most cases here near identical.']);
disp(['This demo will generate a new set of data each time it is run. Occassionally,']);
disp(['fits may fail.']);

%Nelder-Mead search options
options = PAL_minimize('options');  %decrease tolerance (i.e., increase
options.TolX = 1e-09;              %precision).
options.TolFun = 1e-09;
options.MaxIter = 10000;
options.MaxFunEvals = 10000;

%prepare plot
figure('units','pixels','position',[100 100 400 200]);
h = axes('units','pixels','position',[50 50 300 100]);
set(h,'xtick',[-10 -2:1:2 10],'xlim',[-11 11],'ylim',[0 1],'xticklabel',{'-10','-2','-1','0','1','2','10'});
xlabel('Stimulus Intensity');
ylabel('prop correct');
hold on;

%Stimulus intensities. Generating logistic (F) evaluates to near unity
%at 10
StimLevels = [-10 -2 -1 0 1 2 10];

OutOfNum = [50 50 50 50 50 50 50];    %N need not be equal

PF = @PAL_Logistic;             %PF function
paramsValues = [0 1 .05 .05];   %generating values

%Simulate observer
NumPos = PAL_PF_SimulateObserverParametric(paramsValues,StimLevels,OutOfNum,PF);
plot(h,StimLevels,NumPos./OutOfNum,'ko','markersize',6,'markerfacecolor','k');

%Fit PF

searchGrid.alpha = [-1:.05:1];      %structure defining grid to
searchGrid.beta = 10.^[-1:.05:2];   %search for initial values
searchGrid.gamma = [0:.005:.1];     %type help PAL_PFML_Fit for more information
searchGrid.lambda = [0:.005:.1];

paramsFree = [1 1 1 1]; %[threshold slope guess lapse] 1: free, 0:fixed
                        %Note:
                        %Entry for paramsFree(3) (i.e., guess rate) will be
                        %ignored when guess and lapse rates are contrained
                        %to be equal below.

%Fit PF, allowing guess and lapse rate to take on different values
[paramsFitted LL exitflag] = PAL_PFML_Fit(StimLevels, NumPos, OutOfNum, searchGrid, paramsFree, PF,...
    'lapseLimits',[0 1],'searchOptions',options);

if ~exitflag
    disp('Psychometric Function Fit failed! Exiting...');
    break
end

%Add fitted model to plot
plot(h,-10.5:.01:10.5,PF(paramsFitted,-10.5:.01:10.5),'-','color',[0 .7 0],'linewidth',1);

%report parameter estimates
disp(sprintf('\n'));
disp(['Parameter estimates when fitted with gamma and lambda not']);
disp(['constrained to be equal:']);
disp(['Threshold: ' num2str(paramsFitted(1),'%4.3f')]);
disp(['Slope: ' num2str(paramsFitted(2),'%4.3f')]);
disp(['Guess Rate: ' num2str(paramsFitted(3),'%4.3f')]);
disp(['Lapse Rate: ' num2str(paramsFitted(4),'%4.3f')]);
disp(sprintf('\n'));    

%Fit again constraining guess rate and lapse rate to be equal
gammaEQlambda = logical(true); %or: gammaEQlambda = 1;

%Fit model with 'gammaEQlambda' set to TRUE (guess and lapse rate estimates will be equal:
%(Note: the vector assigned to searchGrid.gamma above will now be ignored).
[paramsFitted LL exitflag] = PAL_PFML_Fit(StimLevels, NumPos, OutOfNum, searchGrid, ...
    paramsFree, PF,'lapseLimits',[0 1],'searchOptions',options,'gammaEQlambda',gammaEQlambda);

if ~exitflag
    disp('Psychometric Function Fit failed! Exiting...');
    break
end

%Add fitted model to plot
plot(h,-10.5:.01:10.5,PF(paramsFitted,-10.5:.01:10.5),'-','color',[.7 0 0],'linewidth',1);

%report parameter estimates
disp(sprintf('\n'));
disp(['Parameter estimates when fitted with gamma and lambda']);
disp(['constrained to be equal:']);
disp(['Threshold: ' num2str(paramsFitted(1),'%4.3f')]);
disp(['Slope: ' num2str(paramsFitted(2),'%4.3f')]);
disp(['Guess Rate: ' num2str(paramsFitted(3),'%4.3f')]);
disp(['Lapse Rate: ' num2str(paramsFitted(4),'%4.3f')]);
disp(sprintf('\n'));    

%Fit model with 'gammaEQlambda' set to TRUE and 'lapseFit' set to 'iAPLE':
%'iAPLE': Highest stimulus intensities (and lowest when gammaEQlambda is set
%to TRUE are assumed to be at asymptotic levels. Negative responses at the 
%highest intensity (and positive responses at the lowest intensity) are
%assumed to have happened due to lapses exclusively. Lapses are fit based
%on highest (and lowest) intensities, then lapse (and guess) rate is fixed
%while threshold and slope are estimated based on responses at other
%intensities.
%(the vector assigned to searchGrid.gamma above will now be ignored).
[paramsFitted LL exitflag] = PAL_PFML_Fit(StimLevels, NumPos, OutOfNum, searchGrid, ...
    paramsFree, PF,'lapseLimits',[0 1],'searchOptions',options,'gammaEQlambda',gammaEQlambda,'lapseFit','iAPLE');

if ~exitflag
    disp('Psychometric Function Fit failed! Exiting...');
    break
end

%Add fitted model to plot
plot(h,-10.5:.01:10.5,PF(paramsFitted,-10.5:.01:10.5),'-','color',[0 0 0.7],'linewidth',1);

%report parameter estimates
disp(sprintf('\n'));
disp(['Parameter estimates when fitted with gamma and lambda constrained']);
disp(['to be equal and ''lapseFit'' set to ''iAPLE'':']);
disp(['Threshold: ' num2str(paramsFitted(1),'%4.3f')]);
disp(['Slope: ' num2str(paramsFitted(2),'%4.3f')]);
disp(['Guess Rate: ' num2str(paramsFitted(3),'%4.3f')]);
disp(['Lapse Rate: ' num2str(paramsFitted(4),'%4.3f')]);
disp(sprintf('\n'));    

question = sprintf('\nDo you wish to perform a non-parametric bootstrap to obtain \nSEs for the parameter estimates [y/n]? ');
wish = input(question,'s');
if strcmpi(wish,'y')
    B = input('Type number of desired simulations B (try low, things take a while): ');

    [SD paramsSim LLSim converged] = PAL_PFML_BootstrapNonParametric(StimLevels, NumPos, OutOfNum, ...
    [],paramsFree, B, PF, 'lapseLimits',[0 1], 'searchGrid', searchGrid,'gammaEQlambda',gammaEQlambda,'lapseFit','iAPLE');

    disp([ sprintf('\n') 'SE for the threshold: ' num2str(SD(1),'%4.3f')]);
    disp(['SE for the slope: ' num2str(SD(2),'%4.3f')]);
    disp(['SE for the guess rate as well as lapse rate: ' num2str(SD(4),'%4.4f') sprintf('\n')]);
end

question = sprintf('\nDo you wish to perform a parametric bootstrap to obtain \nSEs for the parameter estimates [y/n]? ');
wish = input(question,'s');
if strcmpi(wish,'y')
    B = input('Type number of desired simulations B (try low, things take a while): ');

    [SD paramsSim LLSim converged] = PAL_PFML_BootstrapParametric(StimLevels, OutOfNum, ...
    paramsFitted, paramsFree, B, PF, 'lapseLimits',[0 1], 'searchGrid', searchGrid,'gammaEQlambda',gammaEQlambda,'lapseFit','iAPLE');

    disp([ sprintf('\n') 'SE for the threshold: ' num2str(SD(1),'%4.3f')]);
    disp(['SE for the slope: ' num2str(SD(2),'%4.3f')]);
    disp(['SE for the guess rate as well as the lapse rate: ' num2str(SD(4),'%4.4f') sprintf('\n')]);
end

question = sprintf('\nDo you wish to determine the Goodness-of-Fit for the fit [y/n]? ');
wish = input(question,'s');
if strcmpi(wish,'y')
    B = input('Type number of desired simulations B (try low, things take a while): ');

    [Dev pDev DevSim converged] = PAL_PFML_GoodnessOfFit(StimLevels, NumPos, OutOfNum, ...
    paramsFitted, paramsFree, B, PF, 'searchGrid', searchGrid, 'lapseLimits',[0 1],'gammaEQlambda',gammaEQlambda,'lapseFit','iAPLE');

    disp([sprintf('\n') 'Goodness-of-fit by Monte Carlo: ' num2str(pDev,'%4.4f')])

    if exist('chi2pdf.m') == 2
        disp(['Goodness-of-fit by chi-square approximation: ' num2str(1-chi2cdf(Dev,4),'%4.4f')])
    end
end
