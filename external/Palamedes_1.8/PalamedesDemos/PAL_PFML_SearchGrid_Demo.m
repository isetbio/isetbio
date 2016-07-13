%PAL_PFML_SearchGrid_Demo  Demonstrates Palamedes' brute force
%search for initial PF parameter values to be used in Simplex search.
%
%Note: The description below applies when the demo is run in Matlab. If the
%demo is run under Octave, the demo produces no figures, but instead prints
%results to the command window.
%
%This demo performs two Psychometric Function fits. 
%
%In the first fit, threshold and slope are free parameters, guess and lapse 
%rate are fixed. A search grid containing all 676 combinations of 26 
%discrete threshold values and 26 discrete slope values is constructed. 
%PAL_PFML_Fit first finds the best-fitting PF of the 676 PFs in the search 
%grid, then performs a Simplex search using this best-fitting PF in the grid 
%as the starting point for the Simplex search. A figure is produced which 
%displays the Likelihoods for the 676 PFs in the grid. Also indicated 
%are the location of the best fitting PF of the 676 in the grid as well as 
%that of the best-fitting PF found by the Simplex method. Standard errors
%of the parameter estimates are then found using a parametric as well as a
%non-parametric bootstrap, and the goodness-of-fit is determined. For each 
%of the simulated datasets the best-fitting PF is first found using a 
%brute-force search. This PF is then the starting point for the Simplex 
%search.
%
%In the second fit, threshold, slope and lapse rate are free parameters and
%the guess rate is fixed. A search grid containing all 20,280 combinations 
%of 26 discrete threshold values, 26 discrete slope values, and 30 discrete
%lapse rate values is constructed. PAL_PFML_Fit first finds the best-
%fitting PF of the 20,280 in the search grid, then performs a Simplex 
%search using the best-fitting PF in the grid as the starting point for the 
%Simplex search. A figure is produced which displays the Likelihood values 
%for some of the PFs in the 3-D grid. Specifically: for those PFs that are 
%on the three planes which have constant threshold, slope and lapse rate,
%respectively and intersect at the point in the grid corresponding to the
%best-fitting PF in the grid. Also shown are these 2-D planes individually, 
%each with the location of the best-fitting PF that is in the grid as well 
%as that of the best-fitting PF found by the Simplex method indicated.
%
%Demonstrates usage of Palamedes functions:
%-PAL_PFML_Fit
%-PAL_PFML_BruteForceFit
%-PAL_PFML_BootstrapParametric
%-PAL_PFML_BootstrapNonParametric
%-PAL_PFML_GoodnessOfFit
%secondary:
%PAL_Gumbel
%PAL_PF_SimulateObserverParametric
%
%More information on any of these functions may be found by typing
%help followed by the name of the function. e.g., help PAL_PFML_Fit
%
%NP (March 2011)

clear all

rand('seed', sum(100*clock)); %will produce different results eacht time it is run.

MorO = ~exist('OCTAVE_VERSION');

if ~MorO
    fprintf('\nUnder Octave, Figures do not render (at all!) as intended. Visit\n');
    fprintf('www.palamedestoolbox.org/demosfiguregallery.html to see figures\n');
    fprintf('as intended.\n\n');
end


%First Fit (threshold and slope are free, guess rate and lapse rate are fixed)

NumTrials = 100;
GridGrain = 26; %26 gives nice round numbers. Things will work fine if 
                %changed unless available memory is exceeded.
B = 400;    %Number of simulations to perform for standard error and 
            %Goodness-Of-Fit determinations

PF = @PAL_Gumbel;           %generating and fitted function
paramsGen = [0 1 0.5 0.03]; %generating parameters
StimLevels = PF([paramsGen(1:2) 0 0],[.2 .4 .5 .6 .8 .9 .98],'inverse');
OutOfNum = NumTrials*ones(size(StimLevels));

%Simulate observer
NumPos = PAL_PF_SimulateObserverParametric(paramsGen,StimLevels,OutOfNum,PF);

%Threshold and slope are free parameters, guess and lapse rate are fixed
paramsFree = [1 1 0 0];

%Specify parameter grid to be searched by brute force (2D here since two 
%parameters are free, but may be 0, 1, 2, 3, or 4D)
searchGrid.alpha = -.25:(.25- -.25)/(GridGrain-1):.25; %Values to be searched in brute force grid search
searchGrid.beta = logspace(-.5,.5,GridGrain); %Values to be searched in brute force grid search
searchGrid.gamma = 0.5; %scalar: parameter is fixed, so nothing to search
searchGrid.lambda = .03; %scalar: parameter is fixed, so nothing to search

%Fit function. PAL_PFML_Fit will first perform brute force search through
%grid specified in 'searchGrid', then use result as starting point in Simplex
%search.
paramsValues = PAL_PFML_Fit(StimLevels, NumPos, OutOfNum, searchGrid, paramsFree, PF);

%Visualize things:
%Find out what initial search values were used in Simplex search and
%generate Likelihood values for all points in 'searchGrid'
[paramsBruteForce LL LLspace] = PAL_PFML_BruteForceFit(StimLevels, NumPos, OutOfNum, searchGrid, PF);

if MorO     %running under Matlab

    %Visualize search grid.
    figure('units','pixels','position',[100 100 800 600],'name','Free: Threshold, Slope. Fixed: Guess, Lapse');
    LLspace = PAL_Scale0to1(LLspace)*64;
    image(squeeze(flipud(LLspace(:,:)')))
    set(gca,'units','pixels','position',[75 75 400 400]);
    set(gca,'xtick',[(GridGrain)/10+.5:GridGrain/5:(GridGrain)*.9+.5],'xticklabel',{'-0.2','-0.1','0','0.1','0.2'},'ytick',[(GridGrain)/10+.5:GridGrain/5:(GridGrain)*.9+.5],'yticklabel',{'0.4','0.2','0','-0.2','-0.4'},'fontsize',12);
    xlabel('Threshold');
    ylabel('Log10(Slope)');

    set(gca,'handlevisibility','off');
    plot(paramsBruteForce(1),log10(paramsBruteForce(2)),'ko','markersize',12,'linewidth',2);
    set(gca,'units','pixels','position',[75 75 400 400]);
    axis([-.25 .25 -.5 .5]);
    text(0,.53,'Likelihood search grid','fontsize',15,'horizontalalignment','center')
    axis off;
    hold on
    plot(paramsValues(1),log10(paramsValues(2)),'ko','markerfacecolor','k');

    set(gca,'handlevisibility','off');
    plot(.15,.4,'ko','markersize',12,'linewidth',2);
    hold on
    plot(.15,.15,'ko','markerfacecolor','k');
    set(gca,'units','pixels','position',[25 500 700 100]);
    axis([0 1 0 1]);
    axis off;
    text(.75,.9,'\alpha','fontsize',12,'horizontalalignment','center')
    text(.85,.9,'log10(\beta)','fontsize',12,'horizontalalignment','center')
    message = ['Generating parameters:'];
    text(.7,.65,message,'fontsize',12,'horizontalalignment','right')
    text(.75,.65, '0','horizontalalignment','center','fontsize',12);
    text(.85,.65, '0','horizontalalignment','center','fontsize',12);
    message = ['Initial values for parameters from brute force search:'];
    text(.7,.4,message,'fontsize',12,'horizontalalignment','right')
    text(.75,.4, num2str(paramsBruteForce(1), '%6.4f'),'horizontalalignment','center','fontsize',12);
    text(.85,.4, num2str(log10(paramsBruteForce(2)), '%6.4f'),'horizontalalignment','center','fontsize',12);
    message = ['Final estimates for parameters after simplex search:'];
    text(.7,.15,message,'fontsize',12,'horizontalalignment','right')
    text(.75,.15, num2str(paramsValues(1), '%6.4f'),'horizontalalignment','center','fontsize',12);
    text(.85,.15, num2str(log10(paramsValues(2)), '%6.4f'),'horizontalalignment','center','fontsize',12);


    set(gca,'handlevisibility','off');
    minX = PAL_Gumbel([paramsGen(1:2) 0 0],.1,'inverse');
    maxX = PAL_Gumbel([paramsGen(1:2) 0 0],.999,'inverse');
    plot(StimLevels,NumPos./OutOfNum,'kd','markersize',8,'markerfacecolor','k','linewidth',2);
    hold on
    plot([minX:(maxX-minX)/100:maxX],PF([paramsValues(1) paramsValues(2) .5 .03], [minX:(maxX-minX)/100:maxX]),'k-');
    set(gca,'units','pixels','position',[550 325 200 150]);
    set(gca,'xtick',[-.5:.25:.5])
    set(gca,'ytick',[.5:.1:1])
    xlabel('Stimulus Intensity')
    ylabel('proportion correct')
    axis([PAL_Gumbel([paramsGen(1:2) 0 0],.15,'inverse') PAL_Gumbel([paramsGen(1:2) 0 0],.99,'inverse') 0.4 1]);
    drawnow

    set(gca,'handlevisibility','off');
    plot(-1,-1);
    axis off
    axis([0 1 0 1]);
    set(gca,'units','pixels','position',[500 75 275 200]);
    text(.1,.9,'Standard Errors:','Fontsize',14);
    text(.6,.7,'\alpha','fontsize',12,'horizontalalignment','center');
    text(.8,.7,'log10(\beta)','fontsize',12,'horizontalalignment','center');
    text(.1,.4,'Parametric:');          %corrected (switched) labels in Palamedes v. 1.3.1
    text(.1,.55,'Non-Parametric:');
    text(.1,.25,'Goodness-of-Fit (PF model versus ');
    text(.1,.15,'saturated model p-value):');
    rectangle('position',[.5 .5 .4 .1],'facecolor',[1 1 1],'edgecolor',[1 1 1]);
    rectangle('position',[.5 .35 .4 .1],'facecolor',[1 1 1],'edgecolor',[1 1 1]);
    rectangle('position',[.65 .1 .25 .1],'facecolor',[1 1 1],'edgecolor',[1 1 1]);
    drawnow
    text(.7,.55,'in progress','horizontalalignment','center');
    drawnow

else    %running under Octave
    
    fprintf('\nThreshold and Slope free, Guess and Lapse fixed:\n');
    fprintf('Brute-Force search: \nThreshold: %8.3f, log10(Slope): %8.3f\n', paramsBruteForce(1), log10(paramsBruteForce(2)));
    fprintf('Simplex search: \nThreshold: %8.3f, log10(Slope): %8.3f\n\n', paramsValues(1), log10(paramsValues(2)));
    fprintf('Determining standard errors, please wait.\n');
    
end
    
%Perform Bootstrap. Data from each simulation will be fit with simplex
%after finding initial search values using brute force search. Display
%results as they become available.

[SD paramsSim] = PAL_PFML_BootstrapNonParametric(StimLevels, NumPos, OutOfNum, [], paramsFree, B, PF, 'searchGrid', searchGrid);
SDlog10slope = std(log10(paramsSim(:,2)));

if MorO

    rectangle('position',[.5 .5 .4 .1],'facecolor',[1 1 1],'edgecolor',[1 1 1]);
    text(.6,.55,num2str(SD(1),'%6.3f'),'horizontalalignment','center');
    text(.8,.55,num2str(SDlog10slope,'%6.3f'),'horizontalalignment','center');
    text(.7,.4,'in progress','horizontalalignment','center');
    drawnow
    
else
    
    fprintf('Standard errors by non-parametric bootstrap:\n');
    fprintf('Threshold: %8.3f, log10(Slope): %8.3f\n\n', SD(1), SDlog10slope);
    
end

[SD paramsSim] = PAL_PFML_BootstrapParametric(StimLevels, OutOfNum, paramsValues, paramsFree, B, PF, 'searchGrid', searchGrid);
SDlog10slope = std(log10(paramsSim(:,2)));

if MorO
    
    rectangle('position',[.5 .35 .4 .1],'facecolor',[1 1 1],'edgecolor',[1 1 1]);
    text(.6,.4,num2str(SD(1),'%6.3f'),'horizontalalignment','center');
    text(.8,.4,num2str(SDlog10slope,'%6.3f'),'horizontalalignment','center');
    text(.775,.15,'in progress','horizontalalignment','center');
    drawnow

else
    
    fprintf('Standard errors by parametric bootstrap:\n');
    fprintf('Threshold: %8.3f, log10(Slope): %8.3f\n\n', SD(1), SDlog10slope);
    fprintf('Determining Goodness-of-fit, please wait.\n');
    
end    
    
%Perform Goodness-of-Fit

[Dev pDev DevSim converged] = PAL_PFML_GoodnessOfFit(StimLevels, NumPos, OutOfNum, paramsValues, paramsFree, B, PF, 'searchGrid', searchGrid);

if MorO

    rectangle('position',[.65 .1 .25 .1],'facecolor',[1 1 1],'edgecolor',[1 1 1]);
    text(.775,.15,num2str(pDev,'%6.3f'),'horizontalalignment','center');

else
    
    fprintf('Goodness-of-fit:\n');
    fprintf('p = %8.3f\n\n', pDev);
    
end

%Second Fit (threshold, slope, lapse rate are free, guess rate fixed)

NumTrials = 1000; %Lots of trials to get meaningful lapse rate estimate
GridGrain = 26; %26 gives nice round numbers
GridGrainLambda = 30; %30 gives nice round numbers

PF = @PAL_Gumbel;
paramsGen = [0 1 0.5 0.03];
StimLevels = PF([paramsGen(1:2) 0 0],[.2 .4 .5 .6 .8 .9 .98],'inverse');
OutOfNum = NumTrials*ones(size(StimLevels));
NumPos = PAL_PF_SimulateObserverParametric(paramsGen,StimLevels,OutOfNum,PF);

paramsFree = [1 1 0 1];
lapseLimits = [0 0.06];

minAlpha = -.25;
maxAlpha = .25;

searchGrid.alpha = minAlpha:(maxAlpha-minAlpha)/(GridGrain-1):maxAlpha;
searchGrid.beta = logspace(-.5,.5,GridGrain);
searchGrid.gamma = 0.5;
searchGrid.lambda = 0.002:(.06-.002)/(GridGrainLambda-1):.06; %0 left out here only to make visualization work, otherwise 0 may/should be included

paramsValues = PAL_PFML_Fit(StimLevels, NumPos, OutOfNum, searchGrid, paramsFree, PF,'lapselimits',lapseLimits);

%Visualize things
[paramsBruteForce LL LLspace] = PAL_PFML_BruteForceFit(StimLevels, NumPos, OutOfNum, searchGrid, PF);

if MorO

    figure('units','pixels','position',[100 100 800 600],'name','Free: Threshold, Slope, Lapse. Fixed: Guess');
    [paramsGrid.alpha paramsGrid.beta paramsGrid.gamma paramsGrid.lambda] = ndgrid(searchGrid.alpha,searchGrid.beta,searchGrid.gamma,searchGrid.lambda);

    paramsGrid.alpha = squeeze(paramsGrid.alpha);
    paramsGrid.beta = squeeze(paramsGrid.beta);
    paramsGrid.lambda = squeeze(paramsGrid.lambda);
    LLspace = squeeze(LLspace);

    %set all values other than on three planes to minimum value on the three planes, then use
    %PAL_Scale0to1
    LLspace(paramsGrid.alpha ~= paramsBruteForce(1) & paramsGrid.beta ~= paramsBruteForce(2) & paramsGrid.lambda ~= paramsBruteForce(4)) = min(min(min(LLspace(paramsGrid.alpha == paramsBruteForce(1) | paramsGrid.beta == paramsBruteForce(2) | paramsGrid.lambda == paramsBruteForce(4)))));
    LLspace = PAL_Scale0to1(LLspace);

    paramsGrid.alpha= permute(paramsGrid.alpha,[3 1 2]);
    paramsGrid.beta = permute(paramsGrid.beta,[3 1 2]);
    paramsGrid.lambda = permute(paramsGrid.lambda,[3 1 2]);
    LLspace = permute(LLspace,[3 1 2]);

    slice(paramsGrid.alpha, paramsGrid.lambda, log10(paramsGrid.beta), LLspace,paramsBruteForce(1),paramsBruteForce(4),log10(paramsBruteForce(2)))
    colormap('jet')
    axis([min(searchGrid.alpha) max(searchGrid.alpha) min(searchGrid.lambda) max(searchGrid.lambda) min(log10(searchGrid.beta)) max(log10(searchGrid.beta))])
    xlabel('Threshold (\alpha)','fontsize',12)
    ylabel('Lapse (\lambda)','fontsize',12)
    zlabel('log10(Slope (\beta))','fontsize',12);
    set(gca,'ztick', [-.2:.1:.2])
    set(gca,'ytick', [.002 .01:.01:.06])
    set(gca,'ztick', [-.4:.2:.4])
    set(gca,'units','pixels','position',[75 75 400 340])
    set(gca,'handlevisibility','off');

    %show three planes with brute-force fit and simplex fit of parameter values
    %indicated
    image(flipud(squeeze(LLspace(searchGrid.lambda == paramsBruteForce(4),:,:))')*64)
    set(gca,'units','pixels','position',[600 440 150 150])
    set(gca,'xtick',[[-.2:.1:.2].*(GridGrain-1)/.5+(GridGrain+1)/2],'xticklabel',{'-0.2','-0.1','0','0.1','0.2'},'ytick',[[-.4:.2:.4]*(GridGrain-1)/1+(GridGrain+1)/2],'yticklabel',{'0.4','0.2','0','-0.2','-0.4'},'fontsize',10);
    xlabel('Threshold');
    ylabel('Log10(Slope)');
    set(gca,'handlevisibility','off');
    plot(paramsBruteForce(1),log10(paramsBruteForce(2)),'ko','markersize',12,'linewidth',2);
    set(gca,'units','pixels','position',[600 440 150 150]);
    axis([searchGrid.alpha(1) searchGrid.alpha(GridGrain) log10(searchGrid.beta(1)) log10(searchGrid.beta(GridGrain))]);
    axis off;
    hold on
    plot(paramsValues(1),log10(paramsValues(2)),'ko','markerfacecolor','k');
    set(gca,'handlevisibility','off');

    image(flipud(squeeze(LLspace(:,:,searchGrid.beta == paramsBruteForce(2))))*64)
    set(gca,'units','pixels','position',[600 245 150 150])
    set(gca,'xtick',[[-.2:.1:.2].*(GridGrain-1)/.5+(GridGrain+1)/2],'xticklabel',{'-0.2','-0.1','0','0.1','0.2'},'ytick',[[.062-([.06:-.01:.01 .002])]*(GridGrainLambda-1)/(.06-.002)],'yticklabel',{'0.06','0.05','0.04','0.03','0.02','0.01','0.002'},'fontsize',10);
    xlabel('Threshold');
    ylabel('Lapse');
    set(gca,'handlevisibility','off');
    plot(paramsBruteForce(1),paramsBruteForce(4),'ko','markersize',12,'linewidth',2);
    set(gca,'units','pixels','position',[600 245 150 150]);
    axis([searchGrid.alpha(1) searchGrid.alpha(GridGrain) searchGrid.lambda(1) searchGrid.lambda(GridGrainLambda)]);
    axis off;
    hold on
    plot(paramsValues(1),paramsValues(4),'ko','markerfacecolor','k');
    set(gca,'handlevisibility','off');

    image(fliplr(flipud(squeeze(LLspace(:,searchGrid.alpha == paramsBruteForce(1),:))'))*64)
    set(gca,'units','pixels','position',[600 50 150 150])
    set(gca,'xtick',[[.062-([.06:-.02:.02 .002])]*(GridGrainLambda-1)/(.06-.002)],'xticklabel',{'0.06','0.04','0.02','0.002'},'ytick',[[-.4:.2:.4]*(GridGrain-1)/1+(GridGrainLambda+1)/2],'yticklabel',{'0.4','0.2','0','-0.2','-0.4'},'fontsize',10);
    xlabel('Lapse');
    ylabel('log10(Slope)');
    set(gca,'handlevisibility','off');
    plot(.06-paramsBruteForce(4),log10(paramsBruteForce(2)),'ko','markersize',12,'linewidth',2);
    set(gca,'units','pixels','position',[600 50 150 150]);
    axis([searchGrid.lambda(1) searchGrid.lambda(GridGrainLambda) log10(searchGrid.beta(1)) log10(searchGrid.beta(GridGrain))]);
    %text(0,.53,'Likelihood search grid','fontsize',15,'horizontalalignment','center')
    axis off;
    hold on
    plot(.06-paramsValues(4),log10(paramsValues(2)),'ko','markerfacecolor','k');
    set(gca,'handlevisibility','off');

    %Show data points and fitted PF
    minX = PAL_Gumbel([paramsGen(1:2) 0 0],.1,'inverse');
    maxX = PAL_Gumbel([paramsGen(1:2) 0 0],.999,'inverse');
    plot(StimLevels,NumPos./OutOfNum,'kd','markersize',8,'markerfacecolor','k','linewidth',2);
    hold on
    plot([minX:(maxX-minX)/100:maxX],PF([paramsValues(1) paramsValues(2) .5 .03], [minX:(maxX-minX)/100:maxX]),'k-');
    set(gca,'units','pixels','position',[50 435 200 150]);
    set(gca,'xtick',[-.5:.25:.5])
    set(gca,'ytick',[.5:.1:1])
    xlabel('Stimulus Intensity')
    ylabel('proportion correct')
    axis([PAL_Gumbel([paramsGen(1:2) 0 0],.15,'inverse') PAL_Gumbel([paramsGen(1:2) 0 0],.99,'inverse') 0.4 1]);
    set(gca,'handlevisibility','off');

    %Display fitted parameter values
    plot(.35,.55,'ko','markersize',12,'linewidth',2);
    set(gca,'units','pixels','position',[260 435 280 150]);
    axis([0 1 0 1]);
    axis off;
    text(.02,.9,'Parameters:','Fontsize',14)
    text(.02,.7,'Generating:','Fontsize',10)
    text(.02,.55,'Brute-Force:')
    text(.5,.9,'\alpha','fontsize',12,'horizontalalignment','center')
    text(.7,.9,'log10(\beta)','fontsize',12,'horizontalalignment','center')
    text(.9,.9,'\lambda','fontsize',12,'horizontalalignment','center')
    rectangle('position',[.4 0 .6 .8],'facecolor',[1 1 1],'edgecolor',[1 1 1]);
    text(.5,.7,'0','fontsize',10,'horizontalalignment','center')
    text(.7,.7,'0','fontsize',10,'horizontalalignment','center')
    text(.9,.7,'0.03','fontsize',10,'horizontalalignment','center')
    text(.5,.55,num2str(paramsBruteForce(1),'%6.4f'),'fontsize',10,'horizontalalignment','center')
    text(.7,.55,num2str(log10(paramsBruteForce(2)),'%6.4f'),'fontsize',10,'horizontalalignment','center')
    text(.9,.55,num2str(paramsBruteForce(4),'%6.4f'),'fontsize',10,'horizontalalignment','center')
    hold on
    plot(.35,.4,'ko','markerfacecolor','k');
    text(.02,.4,'Simplex:')
    text(.5,.4,num2str(paramsValues(1),'%6.4f'),'fontsize',10,'horizontalalignment','center')
    text(.7,.4,num2str(log10(paramsValues(2)),'%6.4f'),'fontsize',10,'horizontalalignment','center')
    text(.9,.4,num2str(paramsValues(4),'%6.4f'),'fontsize',10,'horizontalalignment','center')
    text(.02,.25,'SE (parametric):')
    text(.7,.25,'in progress','horizontalalignment','center')
    text(.02,.1,'Goodness-Of-Fit:')
    drawnow

else

    fprintf('\n\nThreshold, Slope, and Lapse free, Guess fixed:\n');
    fprintf('Brute-Force search: \nThreshold: %8.3f, log10(Slope): %8.3f, Lapse: %8.3f\n', paramsBruteForce(1), log10(paramsBruteForce(2)),paramsBruteForce(4));
    fprintf('Simplex search: \nThreshold: %8.3f, log10(Slope): %8.3f, Lapse: %8.3f\n\n', paramsValues(1), log10(paramsValues(2)), paramsValues(4));
    fprintf('Determining standard errors, please wait.\n');
    
end

[SD paramsSim] = PAL_PFML_BootstrapParametric(StimLevels, OutOfNum, paramsValues, paramsFree, B, PF, 'searchGrid', searchGrid,'lapseLimits',lapseLimits);

if MorO

    rectangle('position',[.4 .15 .6 .15],'facecolor',[1 1 1],'edgecolor',[1 1 1]);
    text(.5,.25,num2str(SD(1),'%6.4f'),'horizontalalignment','center');
    SDlog10slope = std(log10(paramsSim(:,2)));
    text(.7,.25,num2str(SDlog10slope,'%6.4f'),'horizontalalignment','center');
    text(.9,.25,num2str(SD(4),'%6.4f'),'horizontalalignment','center');
    text(.7,.1,'in progress','horizontalalignment','center')
    drawnow

else

    fprintf('Standard errors by parametric bootstrap:\n');
    fprintf('Threshold: %8.3f, log10(Slope): %8.3f, Lapse: %8.3f\n\n', SD(1), SDlog10slope,SD(4));    
    fprintf('Determining Goodness-of-fit, please wait.\n');
    
end

[Dev pDev DevSim converged] = PAL_PFML_GoodnessOfFit(StimLevels, NumPos, OutOfNum, paramsValues, paramsFree, B, PF, 'searchGrid', searchGrid,'lapseLimits',lapseLimits);

if MorO

    rectangle('position',[.4 0 .6 .15],'facecolor',[1 1 1],'edgecolor',[1 1 1]);
    text(.5,.1,num2str(pDev,'%6.4f'),'horizontalalignment','center');

else
    
    fprintf('Goodness-of-fit:\n');
    fprintf('p = %8.3f\n\n', pDev);

end