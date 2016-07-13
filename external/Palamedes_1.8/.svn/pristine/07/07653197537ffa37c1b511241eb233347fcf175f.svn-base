%
%PAL_PFLR_CustomDefine_Demo Demonstrates custom-defined reparametrization 
%of the parameters of a PF across multiple datasets. Custom
%reparametrization was introduced in Palamedes version 1.1.0.
%
%In this example, data are generated for two 'tasks', 8 sessions for each
%task. The generating thresholds vary according to exponential decay 
%functions which are identical in the two tasks (the exponential decay 
%function is defined in the function ParametrizeThresholds.m, which needs 
%to exist in the folder from which this demo is run. The code for this 
%function is given below). Stimulus selection is governed by the Psi Method
%(Kontsevich & Tyler, 1999).
%
%As data are collected, thresholds estimates for each session are
%plotted along with estimated standard errors (estimates based on Bayesian 
%method, with apologies to those offended by our mixing of Maximum 
%Likelihood (ML) and Bayesian methods in this demo. It would have been more 
%proper, but also more time-consuming, to use ML estimates in the plot).
%
%Data are then fit under two models. The fuller model assumes that in all 
%sessions the probabilities of a correct response as a function of stimulus 
%intensity is governed by a Gumbel function. It also assumes that 
%thresholds follow an exponential decay function both in Task A and Task B. 
%It assumes that the slopes for the PFs in all sessions of Task A are 
%equal. It assumes that in all sessions of Task B slopes are equal also.
%It assumes the guess-rate for all PFs equals 0.5, and the lapse rate 
%for all PFs equals 0.02. 
%   Since the exponential decay function is characterized by three 
%parameters, the fuller model is an 8 parameter model: 3 parameters for the 
%exponential decay function under Task A, 3 for the exponential decay 
%function under Task B, 1 slope under Task A, 1 slope under Task B.
%
%The lesser model is similar to the fuller except that under the lesser
%model the exponential decay functions under both tasks are identical and
%the slope of PFs do not differ between Task A and Task B. Thus the lesser
%model has 4 free parameters: 3 parameters for the (single, shared) 
%exponential decay function and 1 shared slope.
%
%PAL_PFLR_Modelcomparison is used to compare the two models statistically.
%Also calculated are AIC (Akiake's Information Criterion) values for both 
%models. Generally, the model comparison should favor the lesser model 
%(after all, it is the model under which the dataset is generated). 
%
%Next, a goodness-of-fit test of the lesser model is performed. The lesser 
%model is compared against the saturated model. Generally, this test
%should result in an acceptable goodness-of-fit (again, the dataset is, in
%fact, generated according to the lesser model).
%
%If the demo is run on a machine on which the Matlab statistics toolbox
%is installed, the theoretical chi-square distributions will be shown along
%the empirical sampling distribution of TLR values. Note that the
%systematic difference between the empirical and theoretical sampling
%distributions is to be expected, as data collection was governed by an
%adaptive method (e.g., Wichmann & Hill, 2001).
%
%Note that besides the function ParametrizeThresholds.m, this demo also
%requires a function ParametrizeSlopes.m. These functions are:
%
%function alphas = ParametrizeThresholds(params)
% 
% alphas(1:8) = (params(1)-params(4)/2) + ...
%   (params(2)-params(5)/2)*exp(-(params(3)-params(6)/2)*[0:7]);
% alphas(9:16) = (params(1)+params(4)/2) + ...
%   (params(2)+params(5)/2)*exp(-(params(3)+params(6)/2)*[0:7]);
%
%
%function betas = ParametrizeSlopes(params)
%
% betas(1:8) = params(1)-params(2);
% betas(9:16) = params(1)+params(2);
%
%
%Demonstrates usage of Palamedes functions:
%-PAL_PFML_FitMultiple
%-PAL_PFML_BootstrapParametricMultiple
%-PAL_PFML_BootstrapNonParametricMultiple
%-PAL_PFML_GoodnessOfFitMultiple
%-PAL_PFLR_ModelComparison
%The above are demonstrated here especially with regard to custom 
%parametrization of parameters as introduced in Palamedes version 1.1.0.
%
%More information on any of these functions may be found by typing
%help followed by the name of the function. e.g., help PAL_PFML_FitMultiple
%
%More information on custom parametrization may be found by typing
%'help PAL_PFML_CustomDefine'
%
%NP (November 2009)

clear all;

if exist('OCTAVE_VERSION');
    fprintf('\nUnder Octave, Figure does not render as intended. Visit\n');
    fprintf('www.palamedestoolbox.org/demosfiguregallery.html to see figure\n');
    fprintf('as intended.\n\n');
    fprintf('Also, this demo takes a long time to complete (especially under\n');
    fprintf('Octave).\n\n');
    cont = input('Do you wish to continue? [y/n]: ', 's');
    if strcmp(cont,'y')
        fprintf('\n\n');
    else
        break;
    end
end

rand('state',sum(100*clock));

message = ['This demo will generate a different set of data each time'...
    ' it is run.'];
disp(message);
message = ['As a result, occassionally, it may result in failed fit(s).'];
disp(message);
message = 'Parametric Bootstrap (1) or Non-Parametric Bootstrap? (2): ';
ParOrNonPar = input(message);
message = ['Try low numbers on next questions first. Things take a while ....'];
disp(message);
message = sprintf('Number of simulations to perform to determine standar');
message = strcat(message, 'd errors: ');
Bse = input(message);
message = sprintf('Number of simulations to perform to determine model c');
message = strcat(message, 'omparison p-values: ');
Bmc = input(message);

numsessions = 8;    %sessions per task
numtrials = 500;    %trials per session

ExpDecayParamsMean = [1 2 .5];  %theta 1 thru 3 in figure
ExpDecayParamsDiff = [0 0 0];   %theta 4 thru 6 in figure

ExpDecayParameters = [ExpDecayParamsMean ExpDecayParamsDiff];

%see text above for code of ParametrizeThresholds function
Thresholds = ParametrizeThresholds(ExpDecayParameters); %Generating
slope = .5;     %Generating value
guess = .5;     %Generating value
lapse = .02;    %Generating value

PF = @PAL_Gumbel;   %Generating PF

options = PAL_minimize('options');   %Nelder-Mead search options
options.MaxFunEvals = 5000;         %Allow more function evals
options.MaxIter = 5000;             %Allow more iterations
options.TolFun = 1e-08;             %Increase desired precision on LL
options.TolX = 1e-08;               %Increase desired precision on params
options.Display = 'off';            %suppress PAL_minimize messages


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate data using Psi Method procedure:

%Define prior
priorAlphaRange = 0:.2:4; %values of alpha to include in prior
priorBetaRange = -1:.05:1;  %values of beta to include in prior

%Stimulus values to select from
stimRange = [0:.5:6];   
                            
%2-D Gaussian prior
[a b] = ndgrid(priorAlphaRange,priorBetaRange);
prior = PAL_pdfNormal(a,1,2).*PAL_pdfNormal(b,0,1);
clear a;
clear b;
prior = prior./sum(sum(prior)); %prior should sum to 1


%To be assumed during Psi Method procedure
gamma = 0.5;            %Guess rate
lambda = .02;           %Lapse rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%PLOTTING CODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

title = ['Model fitting and model comparison using non-linear '...
    'reparametrization'];

    figure('name',title,'units','pixels','position',[100 100 1070 625]);

    axes('units','pixels','position',[50 325 225 175]);
    hold on
    h1 = gca;
    axis(h1,[.5 8.5 0 4]);
    set(h1, 'fontsize',12);
    set(h1, 'Xtick',[1:8]);
    xlabel('Session');
    ylabel('Threshold ( \alpha)');
    text(6,3.5,'Task A','fontsize',12);
    text(.2,6.5,'\bullet Individual fits (Psi Method)');

    axes('units','pixels','position',[325 325 225 175]);
    h2 = gca;
    axis(h2,[.5 8.5 0 4]);
    set(h2, 'fontsize',12);
    set(h2, 'Xtick',[1:8]);
    xlabel('Session');
    ylabel('Threshold ( \alpha)');
    text(6,3.5,'Task B','fontsize',12);

    axes('units','pixels','position',[50 50 500 175]);
    h3 = gca;
    set(h3, 'fontsize',12);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%END PLOTTING CODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for session = 1:16  %8 sessions under 'task A' 8 under 'task B'

    PM = PAL_AMPM_setupPM('priorAlphaRange',priorAlphaRange,...
    'priorBetaRange',priorBetaRange, 'numtrials',numtrials, 'PF' , PF,...
    'prior',prior,'stimRange',stimRange,'gamma',gamma,'lambda',lambda);

    %Trial loop

    while PM.stop ~= 1

        %simulate trial
        
        response = rand(1) < PF([Thresholds(session) slope guess lapse],...
            PM.xCurrent);
        PM = PAL_AMPM_updatePM(PM,response);

        %Show updated posterior
        axes(h3)
        contour(squeeze(PM.pdf)',15);        
        set(gca,'xtick',[1:5:21],'xticklabel',{'0','1','2','3','4'});
        xlabel('Threshold');
        set(gca,'ytick',[1:20:41],'yticklabel',{'-1','0','1'});
        ylabel('log10(Slope)');
        text(1,45,'Psi Method Posterior:','FontSize',12)
        drawnow
        
    end
    NumPos(session,:) = PM.response;
    StimLevels(session,:) = PM.x(1:numtrials);
    ThresholdFit(session) = ...
        PM.thresholdUniformPrior(length(PM.thresholdUniformPrior));
    seThresholdFit(session) = ...
        PM.seThresholdUniformPrior(length(PM.seThresholdUniformPrior));
    if session < 9
        axes(h1)
        x = session;
        hold on
    else
        axes(h2)
        x = session - 8;
        hold on
    end
    plot(x,ThresholdFit(session),'o','MarkerSize',4,'color',[0 0 0],...
        'MarkerFaceColor',[0 0 0]);
    line([x x],[ThresholdFit(session)-seThresholdFit(session) ...
        ThresholdFit(session)+seThresholdFit(session)],'linewidth',2,...
        'color',[0 0 0]);
    drawnow
end    
    
OutOfNum = ones(size(NumPos));

[StimLevels NumPos OutOfNum] = PAL_PFML_GroupTrialsbyX(StimLevels, ...
    NumPos, OutOfNum); %within session, combine trials of equal intensity

numParamsSat = sum(sum(OutOfNum~=0));  %Number of parameters for saturated 
                                       %model

%Set up reparametrization structure for lesser model
funcParamsL = PAL_PFML_setupParametrizationStruct;

%lesser model forces identical learning curves for Task A and Task B
%and forces all 16 PFs to have identical slopes
funcParamsL.funcA = @ParametrizeThresholds; %custom-written function (see
                                        %help comments above for code)

funcParamsL.paramsValuesA = [1 2 .5 0 0 0]; %guesses/fixed values
funcParamsL.paramsFreeA = [1 1 1 0 0 0];    %free/fixed parameters
%fixing theta 4 thru 6 at 0 forces identical learning curves for Task A 
%and Task B

%lesser model forces identical slopes for all 16 PFs (see help comments
%above for code of ParametrizeSlopes).
funcParamsL.funcB = @ParametrizeSlopes;
funcParamsL.paramsValuesB = [.5 0];
funcParamsL.paramsFreeB = [1 0];

params(1:16,1) = 0;     %Will be ignored (thresholds are reparametrized)
params(1:16,2) = .5;    %Will be ignored (slopes are reparametrized)
params(1:16,3) = .5;    %fixed value guess rate
params(1:16,4) = .02;   %fixed value lapse rate

%Fit lesser model
[paramsL LLL trash trash funcParamsL numParamsL] = ...
    PAL_PFML_FitMultiple(StimLevels, NumPos, OutOfNum, params, PF, ...
    'thresholds', funcParamsL,'slopes', funcParamsL,'lapserates','fixed',...
    'searchOptions',options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%PLOTTING CODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    x = 1:.01:8;
    y = (funcParamsL.paramsValuesA(1)-funcParamsL.paramsValuesA(4)) + ...
        (funcParamsL.paramsValuesA(2)-funcParamsL.paramsValuesA(5))*...
        exp(-(funcParamsL.paramsValuesA(3)-funcParamsL.paramsValuesA(6))...
        *(x-1));
    axes(h1)
    plot(x,y,'-','linewidth',2,'color',[.7 0 0]);
    MT = ['PF Threshold: \alpha(m, session) = (\theta_1 + m\times'...
    '\theta_4)+(\theta_2 + m\times\theta_5)\timesexp(-(session-1)'...
    '\times(\theta_3 + m\times\theta_6))'];
    text(0.5,6,MT,'color',[0 0 0]);
    text(0.5,5.5,'PF slope: \beta(m) = \omega_1 + m\times\omega_2',...
        'color',[0 0 0]);
    s = warning('off','MATLAB:gui:latexsup:BadTeXString');
    text(19,5.75,'}','FontSize',35);
    drawnow
    warning(s); %returns warnings to previous state
    text(20,6, 'Task A: m = -1');
    text(20,5.5, 'Task B: m = +1');
    MT = ['Lesser Model: \theta_1, \theta_2, \theta_3, \omega_1 are free'...
        ', but \theta_4 = \theta_5 = \theta_6 = \omega_2 = 0. AIC: ',...
        num2str(-2*LLL+2*4,'%6.3f')];
    text(0.5,5,MT,'color',[0.7 0 0]);
    axes(h2)
    plot(x,y,'-','linewidth',2,'color',[.7 0 0]);
    set(h1,'HandleVisibility','off');
    set(h2,'HandleVisibility','off');
    
    clf;

    set(h1,'HandleVisibility','on');
    set(h2,'HandleVisibility','on');

    axes('units','pixels','position',[50 50 165 175]);
    plot(1.1:3.1, funcParamsL.paramsValuesA(1:3), 'o', 'markersize',4,...
        'markerfacecolor',[.7 0 0],'color',[.7 0 0],'linewidth',2);
    hold on
    h3 = gca;
    axis(h3,[.5 3.5 0 3]);
    set(h3, 'fontsize',12);
    set(h3, 'Xtick',[1:3],'xticklabel',[]);
    text(.5,3.24,'Model Parameters:','FontSize',12);
    text(1,-.25,'\theta_1','HorizontalAlignment','center','FontSize',16);
    text(2,-.25,'\theta_2','HorizontalAlignment','center','FontSize',16);
    text(3,-.25,'\theta_3','HorizontalAlignment','center','FontSize',16);
    ylabel('Value');
    hold on

    axes('units','pixels','position',[245 50 165 175]);
    plot(1.1:3.1, funcParamsL.paramsValuesA(4:6), 'o', 'markersize',4,...
        'markerfacecolor',[.7 0 0], 'color',[.7 0 0],'linewidth',2);
    hold on
    h4 = gca;
    axis(h4,[.5 3.5 -1 1]);
    set(h4, 'fontsize',12);
    set(h4, 'Xtick',[1:3],'xticklabel',[]);
    set(h4, 'Ytick',[-1:1:1]);
    text(1,-7/6,'\theta_4','HorizontalAlignment','center','FontSize',16);
    text(2,-7/6,'\theta_5','HorizontalAlignment','center','FontSize',16);
    text(3,-7/6,'\theta_6','HorizontalAlignment','center','FontSize',16);
    hold on

    axes('units','pixels','position',[440 50 50 175]);
    plot(1.1, funcParamsL.paramsValuesB(1), 'o', 'markersize',4,...
        'markerfacecolor',[.7 0 0], 'color',[.7 0 0],'linewidth',2);
    hold on
    h5 = gca;
    axis(h5,[.5 1.5 .3 .7]);
    set(h5, 'fontsize',12);
    set(h5, 'Xtick',[1:1],'xticklabel',[]);
    set(h5, 'Ytick',[0.3:.2:.7]);
    text(1,.3-1/30,'\omega_1','HorizontalAlignment','center','FontSize',16);
    hold on

    axes('units','pixels','position',[520 50 50 175]);
    plot(1.1, funcParamsL.paramsValuesB(2), 'o', 'markersize',4,...
        'markerfacecolor',[.7 0 0],'color',[.7 0 0],'linewidth',2);
    hold on
    h6 = gca;
    axis(h6,[.5 1.5 -.2 .2]);
    set(h6, 'fontsize',12);
    set(h6, 'Xtick',[1:1],'xticklabel',[]);
    set(h6, 'Ytick',[-.2:.2:.2]);
    text(1,-.2-1/30,'\omega_2','HorizontalAlignment','center','FontSize',16);    
    hold on

    drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%END PLOTTING CODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%Set up reparametrization structure for fuller model
funcParamsF = PAL_PFML_setupParametrizationStruct;

%fuller model allows different curves on thresholds for Task A and Task B
%and allows slopes to differ between Task A and Task B
funcParamsF.funcA = funcParamsL.funcA;
funcParamsF.paramsValuesA = funcParamsL.paramsValuesA;
funcParamsF.paramsFreeA = [1 1 1 1 1 1];
funcParamsF.funcB = funcParamsL.funcB;
funcParamsF.paramsValuesB = funcParamsL.paramsValuesB;
funcParamsF.paramsFreeB = [1 1];

%Fit fuller model
[paramsF LLF trash trash funcParamsF numParamsF] = ...
    PAL_PFML_FitMultiple(StimLevels, NumPos, OutOfNum, paramsL, PF, ...
    'thresholds', funcParamsF,'slopes', funcParamsF,'lapserates','fixed',...
    'searchOptions',options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%PLOTTING CODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    x = 1:.01:8;
    y = (funcParamsF.paramsValuesA(1)-funcParamsF.paramsValuesA(4)) + ...
        (funcParamsF.paramsValuesA(2)-funcParamsF.paramsValuesA(5))*...
        exp(-(funcParamsF.paramsValuesA(3)-funcParamsF.paramsValuesA(6))...
        *(x-1));    
    axes(h1)
    plot(x,y,'--','linewidth',2,'color',[0 .7 0]);
    MT = ['Fuller Model: \theta_1, \theta_2, \theta_3, \theta_4, \theta_5'...
        ', \theta_6, \omega_1, \omega_2 are free. AIC: ' ...
        num2str(-2*LLF+2*8,'%6.3f')];
    text(0.5,4.5,MT,'color',[0 0.7 0]);
    
    x = 1:.01:8;
    y = (funcParamsF.paramsValuesA(1)+funcParamsF.paramsValuesA(4)) + ...
        (funcParamsF.paramsValuesA(2)+funcParamsF.paramsValuesA(5))*...
        exp(-(funcParamsF.paramsValuesA(3)+funcParamsF.paramsValuesA(6))...
        *(x-1));
    axes(h2)
    plot(x,y,'--','linewidth',2,'color',[0 .7 0]);
        
    axes(h3);
    plot(.9:2.9, funcParamsF.paramsValuesA(1:3), 'o','markersize',4,...
        'markerfacecolor',[0 .7 0],'color',[0 .7 0],'linewidth',2);

    axes(h4);
    plot(.9:2.9, funcParamsF.paramsValuesA(4:6), 'o','markersize',4,...
        'markerfacecolor',[0 .7 0], 'color',[0 .7 0],'linewidth',2);

    axes(h5);
    plot(.9, funcParamsF.paramsValuesB(1), 'o','markersize',4,...
        'markerfacecolor',[0 .7 0], 'color',[0 .7 0],'linewidth',2);

    axes(h6);
    plot(.9, funcParamsF.paramsValuesB(2), 'o','markersize',4,...
        'markerfacecolor',[0 .7 0], 'color',[0 .7 0],'linewidth',2);

    drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%END PLOTTING CODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%set up rangeTries structure.
rangeTries.paramsValuesA = [.5 .5 .5 .1 .1 .1];
rangeTries.paramsValuesB = [.5 .1];

maxTries = 20;

%Estimate standard errors for parameters of fuller model using bootstrap
if ParOrNonPar == 1
    [SD paramsSim trash trash SDfunc funcParamsSim] = ...
        PAL_PFML_BootstrapParametricMultiple(StimLevels,OutOfNum, ...
        paramsF, Bse, PF,'thresholds', funcParamsF,'slopes',funcParamsF,...
        'lapserates','fixed', 'maxTries', maxTries, 'rangeTries', ...
        rangeTries,'searchOptions',options);
else
    [SD paramsSim trash trash SDfunc funcParamsSim] = ...
        PAL_PFML_BootstrapNonParametricMultiple(StimLevels,NumPos, ...
        OutOfNum, paramsF, Bse, PF,'thresholds', funcParamsF,'slopes',...
        funcParamsF, 'maxTries', maxTries, 'rangeTries', rangeTries,...
        'searchOptions',options);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%PLOTTING CODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    axes(h3);
    for param = 1:3
        line([param-.1 param-.1], [funcParamsF.paramsValuesA(param)-...
            SDfunc.A(param) funcParamsF.paramsValuesA(param)+...
            SDfunc.A(param)], 'color',[0 .7 0],'linewidth',2);
    end

    axes(h4);
    for param = 4:6
        line([param-3.1 param-3.1], [funcParamsF.paramsValuesA(param)-...
            SDfunc.A(param) funcParamsF.paramsValuesA(param)+...
            SDfunc.A(param)], 'color',[0 .7 0],'linewidth',2);
    end

    axes(h5);
    line([.9 .9], [funcParamsF.paramsValuesB(1)-SDfunc.B(1) ...
        funcParamsF.paramsValuesB(1)+SDfunc.B(1)], 'color',[0 .7 0],...
        'linewidth',2);

    
    axes(h6);
    line([.9 .9], [funcParamsF.paramsValuesB(2)-SDfunc.B(2) ...
        funcParamsF.paramsValuesB(2)+SDfunc.B(2)], 'color',[0 .7 0],...
        'linewidth',2);

    drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%END PLOTTING CODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
%Estimate standard errors for parameters of lesser model using bootstrap
if ParOrNonPar == 1
    [SD paramsSim trash trash SDfunc funcParamsSim] = ...
        PAL_PFML_BootstrapParametricMultiple(StimLevels,OutOfNum, ...
        paramsF, Bse, PF,'thresholds', funcParamsL,'slopes', funcParamsL,...
        'lapserates','fixed', 'maxTries', maxTries, 'rangeTries', ...
        rangeTries,'searchOptions',options);
else
    [SD paramsSim trash trash SDfunc funcParamsSim] = ...
        PAL_PFML_BootstrapNonParametricMultiple(StimLevels,NumPos, ...
        OutOfNum, paramsF, Bse, PF,'thresholds', funcParamsL,'slopes',...
        funcParamsL, 'maxTries', maxTries, 'rangeTries', rangeTries,...
        'searchOptions',options);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%PLOTTING CODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    axes(h3);
    for param = 1:3
        line([param+.1 param+.1], [funcParamsL.paramsValuesA(param)-...
            SDfunc.A(param) funcParamsL.paramsValuesA(param)+...
            SDfunc.A(param)], 'color',[.7 0 0],'linewidth',2);
    end

    axes(h4);
    for param = 4:6
        line([param-2.9 param-2.9], [funcParamsL.paramsValuesA(param)-...
            SDfunc.A(param) funcParamsL.paramsValuesA(param)+...
            SDfunc.A(param)], 'color',[.7 0 0],'linewidth',2);
    end

    axes(h5);
    line([1.1 1.1], [funcParamsL.paramsValuesB(1)-SDfunc.B(1) ...
        funcParamsL.paramsValuesB(1)+SDfunc.B(1)], 'color',[.7 0 0],...
        'linewidth',2);

    axes(h6);
    line([1.1 1.1], [funcParamsL.paramsValuesB(2)-SDfunc.B(2) ...
        funcParamsL.paramsValuesB(2)+SDfunc.B(2)], 'color',[.7 0 0],...
        'linewidth',2);
    
    drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%END PLOTTING CODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%Perform model comparison
[TLR pTLR paramsL paramsF TLRSim converged] = ...
    PAL_PFLR_ModelComparison(StimLevels, NumPos, OutOfNum, paramsL, Bmc,...
    PF,'LesserThresholds',funcParamsL,'FullerThresholds',funcParamsF,...
    'lesserslopes',funcParamsL,'fullerslopes',funcParamsF,'lesserlapse',...
    'fixed','fullerlapse','fixed','maxTries',maxTries,'rangetries',...
    rangeTries,'searchOptions',options);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%PLOTTING CODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Distribution of simulated TLRs
    
    axes('units','pixels','position',[650 325 400 175]);
    [n centers] = hist(TLRSim,40);
    hist(TLRSim,40)
    h = findobj(gca,'Type','patch');
    set(gca,'FOntSize',12)
    set(h,'FaceColor','y','EdgeColor','k')
    set(gca,'xlim',[0 1.2*max(TLR,centers(length(centers)))]);
    xlim = get(gca, 'Xlim');
    hold on
    if exist('chi2pdf.m') == 2
        chi2x = xlim(1):xlim(2)/250:xlim(2);
        [maxim I]= max(n);
        chi2 = chi2pdf(chi2x,numParamsF-numParamsL)*...
            (maxim/chi2pdf(centers(I),numParamsF-numParamsL));
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
        message = ['p_{chi2}: ' num2str(1-chi2cdf(TLR,...
            numParamsF-numParamsL),'%5.4f')];
        text(.95*xlim(2),.7*ylim(2),message,'horizontalalignment','right',...
            'fontsize',10);
    end
    text(0,1.08*ylim(2),'Likelihood Ratio Test: Fuller Model vs. Lesser Model',...
        'color',[0 0 0],'Fontsize',12);
    xlabel('Simulated TLRs','FontSize',12)
    ylabel('frequency','FontSize',12);
    drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%END PLOTTING CODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%Perform goodness-of-fit test of lesser model    
[Dev pDev DevSim converged] = PAL_PFML_GoodnessOfFitMultiple(StimLevels, ...
    NumPos, OutOfNum, paramsL, Bmc, PF,'Thresholds',funcParamsL,'slopes',...
    funcParamsL,'lapserate','fixed','maxtries',maxTries,'rangetries',...
    rangeTries,'searchOptions',options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%PLOTTING CODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Distribution of simulated TLRs
    axes('units','pixels','position',[650 50 400 175]);
    [n centers] = hist(DevSim,40);
    hist(DevSim,40)
    h = findobj(gca,'Type','patch');
    set(gca,'FOntSize',12)
    set(h,'FaceColor','y','EdgeColor','k')
    set(gca,'xlim',[0 1.2*max(Dev,centers(length(centers)))]);
    xlim = get(gca, 'Xlim');
    hold on
    if exist('chi2pdf.m') == 2
        chi2x = xlim(1):xlim(2)/250:xlim(2);
        [maxim I]= max(n);
        chi2 = chi2pdf(chi2x,numParamsSat-numParamsL)*...
            (maxim/chi2pdf(centers(I),numParamsSat-numParamsL));
        plot(chi2x,chi2,'k-','linewidth',2)
    end
    ylim = get(gca, 'Ylim');
    plot(Dev,.05*ylim(2),'kv','MarkerSize',12,'MarkerFaceColor','k')
    text(Dev,.15*ylim(2),'Dev data','Fontsize',11,'horizontalalignment',...
        'center');
    message = ['p_{simul}: ' num2str(pDev,'%5.4f')];
    text(.95*xlim(2),.8*ylim(2),message,'horizontalalignment','right',...
        'fontsize',10);
    if exist('chi2cdf.m') == 2
        message = ['p_{chi2}: ' num2str(1-chi2cdf(Dev,...
            numParamsSat-numParamsL),'%5.4f')];
        text(.95*xlim(2),.7*ylim(2),message,'horizontalalignment','right',...
            'fontsize',10);
    end
    text(0,1.08*ylim(2),'Goodness-of-Fit: Saturated Model vs. Lesser Model',...
        'color',[0 0 0],'Fontsize',12);
    xlabel('Simulated Devs','FontSize',12)
    ylabel('frequency','FontSize',12);
    drawnow
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%END PLOTTING CODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%