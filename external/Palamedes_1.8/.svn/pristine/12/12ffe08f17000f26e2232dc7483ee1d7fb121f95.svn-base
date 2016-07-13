%
%PAL_AMPM_Demo  Demonstrates use of Palamedes routines to implement the
%'psi' adaptive procedure (Kontsevich & Tyler, 1999) and some variations 
%(Prins, 2013 or www.palamedestoolbox.org/psimarginal.html).
%
%Note that many more possibilites exist! (see Prins (2013) or 
%www.palamedestoolbox.org/psimarginal.html)
%
%User is prompted as to whether original psi method, the psi+, or the
%psi-marginal should be demonstrated. The psi-method keeps track of a 2D
%(threshold x slope) posterior and selects stimulus placement that will
%minimize the expected entrop in it. It assumes a fixed lapse and guess 
%rate. The psi+ method is identical except that it keeps track of a 3D 
%posterior that includes the lapse rate. It selects stimulus placement
%such as to reduce expected entropy in the 3D posterior. As a result, many
%trials are spent trying to obtain a precise estimate of the lapse rate.
%The psi-marginal method keeps track of a 3D posterior but selects stimulus
%placement such to reduce expected entropy in the marginal threshold x
%slope posterior (or in the marginal threshold distribution). The former
%optimizes estimation of the threshold and slope parameter while allowing
%the lapse rate to be targeted if it serves the estimation of threshold
%and/or slope. The latter optimizes estimation of threshold only but allows
%slope and lapse rate to be targeted if it helps estimation of the
%threshold parameter. Each of the four versions has its own selection 
%criterion and each will display a distinct stimulus placement strategy 
%that fits that criterion.
%
%When the lapse rate is included in the posterior, the method may have a
%tendency at high N to produce series of consecutive placements at stimulus
%intensities that produce performance near the upper asymptote. This is
%undesirable in that such series are likely to affect an observer's
%vigilance. User will be prompted to indicate whether this should be
%prevented by suspending psi+ or psi-marginal and temporarily act as the
%original psi-method by fixing the lapse rate for a random number of trials
%after a placement at the highest possible intensity. Note that this
%strategy is not fool-proof, especially early in a trial run (the 
%original psi-method produces consecutive series at the extreme amplitudes 
%sometimes as well).
%
%Note that the Psi-marginal method can be set up to treat any of the PF's
%four parameters either as a parameter of primary interest whose
%estimation should be optimized, as a nuisance parameter whose estimation
%should subserve the estimation of the primary parameter, or as a fixed
%parameter.
%
%Note that user may define a prior other than the constrained uniform prior
%used here.
%
%Note that user may constrain the guess rate and the lapse rate to be equal
%as would be appropriate in, say, a vernier acuity task.
%
%Demonstrates usage of Palamedes functions:
%-PAL_AMPM_setupPM
%-PAL_AMPM_updatePM
%secondary:
%PAL_Gumbel
%PAL_findMax
%PAL_PFML_GroupTrialsbyX
%
%More information on any of these functions may be found by typing
%help followed by the name of the function. e.g., help PAL_AMPM_setupPM
%
%References:
%
%Kontsevich, LL & Tyler, CW (1999). Bayesian adaptive estimation of
%psychometric slope and threshold. Vision Research, 39, 2729-2737.
%
%Prins, N (2012a). The psychometric function: The lapse rate revisited.
%Journal of Vision, 12(6):25, 1-16. doi: 10.1167/12.6.25
%
%Prins, N (2012b). The adaptive Psi method and the lapse rate. Poster
%presented at the 12th annual meeting of the Vision Sciences Society.
%http://f1000.com/posters/browse/summary/1090339
%
%Prins, N. (2013). The psi-marginal adaptive method: how to give nuisance 
%parameters the attention they deserve (no more, no less). Journal of
%Vision, 13(7):3, 1-17. doi: 10.1167/13.7.3 
%
%NP (March 2013)

clear all

%opengl software %works without this, but graphs miss axes

if exist('RandStream.m','file');
    s = RandStream.create('mt19937ar','seed','shuffle'); %do something different each time
    RandStream.setGlobalStream(s);
else
    fprintf('\nYour version of Matlab or Octave does not support RandStream:');
    fprintf('\nRandom number generator is not randomly seeded.\n\n');
end
if exist('OCTAVE_VERSION');
    fprintf('\nUnder Octave, Figure will not render exactly (or hardly at all!) as intended. \n');
    fprintf('Visit www.palamedestoolbox.org/demosfiguregallery.html to see figure\n');
    fprintf('as intended.\n\n');
end

marginalize = [];
AvoidConsecutive = 0;

fprintf(1,'\n')
disp('For more information on what goes on in this demo, type help PAL_AMPM_Demo')
disp('or go to: http://www.journalofvision.org/content/13/7/3')
fprintf(1,'\n')
method = input('Original Psi (0), Psi+ (1), or Psi-marginal (2)?: ');
if method == 2
    marginalize = [4];
    if strncmpi(input('Slope of primary interest or nuisance? (p/n): ','s'),'n',1);
        marginalize = [marginalize 2];
    end
end

if method == 1 || method == 2
    fprintf(1,'\n')
    disp('Lengthy consecutive placements at high intensities may be avoided')
    disp('by temporarily using a fixed lapse rate (as in original Psi-method)')
    AvoidConsecutive = strncmpi(input('after a high intensity trial. Do this? (y/n): ','s'),'y',1);
    if AvoidConsecutive
        fprintf(1,'\n')
        disp('After a high intensity trial, the method will assume a fixed')
        disp('lapse rate for a random number of trials. This ''wait time''')
        disp('will be drawn from an exponential mass function in order to')
        disp('maintain constant ''hazard''.')
        fprintf(1,'\n')
        WaitTime = input('Enter average wait time (in number of trials), e.g., 4: ');
    end
end

fprintf(1,'\n')
NumTrials = input('Number of trials: ');
fprintf(1,'\n')

plotGrain = 51; %Allow higher resolution for visualization compared to what Psi-marginal works with
computeGrain = 51; 

suspend = 0;
trackSuspend = [];

gamma = .5; %generating guess rate
lambda = 0.03; %generating lapse rate

PF = @PAL_Gumbel; %generating and assumed function
paramsGen = [0 1 gamma lambda]; 

%Stimulus values the method can select from
stimRange = (linspace(PF([0 1 0 0],.1,'inverse'),PF([0 1 0 0],.9999,'inverse'),21));

%%%%%%%Prepare posterior shown in Figure. This posterior is for
%%%%%%%illustrative purposes only.

alphas = linspace(PF([0 1 0 0],.1,'inverse'),PF([0 1 0 0],.9999,'inverse'),plotGrain);
betas = linspace(log10(.0625),log10(16),plotGrain);
gammas = gamma;
lambdas = [0:.01:.1];
[a, b, l] = ndgrid(alphas,betas,lambdas);

a = permute(a, [2 1 3]);
b = permute(b, [2 1 3]);
l = permute(l, [2 1 3]);

vizLUT = squeeze(PAL_AMPM_CreateLUT(alphas,betas,gammas,lambdas,stimRange,PF,false));
vizLUT = permute(vizLUT,[2 1 3 4]);
posterior = ones(size(a));

%totals keeps track of number of trials presented at each value of stimRange
totals = zeros(size(stimRange));

if ~exist('OCTAVE_VERSION');
    disp('Use space bar and ''g'' to step through or automate, respectively.') 
    kpf = @(src,event) disp('');   %dummy function used in next line
    fh = figure('units','pixels','position',[100 50 1100 700],'keypressfcn',kpf);
else
    fh = figure('units','pixels','position',[100 50 1100 700]);
end

%Something to tag figure labels unto
axes
set(gca,'units','pixels','position',[1 1 1100 700])
set(gca,'xlim',[0 1],'ylim',[0 1])
text(.1,.87,'(a)','fontsize',14)
text(.48,.87,'(b)','fontsize',14)
text(.66,.87,'(c)','fontsize',14)
text(.84,.87,'(d)','fontsize',14)
text(.46,.59,'(e)','fontsize',14)
text(.04,.25,'(f)','fontsize',14)
if method == 0
    text(.03,.97,'(original) Psi method','fontsize',16)
    text(.03,.93,'Figure (a) serves illustrative purposes only. Method maintains posterior in (b) only and selects placement such as to minimize expected entropy in (b).','fontsize',12)
end
if method == 1
    text(.03,.97,'Psi+ method','fontsize',16)
    text(.03,.93,'Method maintains posterior in (a) and selects placement such as to minimize expected entropy in (a).','fontsize',12)
end
if method == 2
    text(.03,.97,'Psi-marginal method','fontsize',16)
    text(.03,.93,'Method maintains posterior in (a) but selects placement such as to minimize expected entropy in (b).','fontsize',12)
end
    
axis off
set(gca,'handlevisibility','off')

%Set up psi(marginal)

alphaOffset = rand(1)*.3 - .15; %jitter prior
betaOffset = rand(1)*.3 - .15;

priorAlphaRange = linspace(PF([0 1 0 0],.1,'inverse'),PF([0 1 0 0],.9999,'inverse'),computeGrain) + alphaOffset;
priorBetaRange =  linspace(log10(.0625),log10(16),computeGrain) + betaOffset;
if method == 0  %psi vs. psi+ and psimarginal
    priorLambdaRange = .03;
else
    priorLambdaRange = [0:.01:.1];
end

%Initialize PM structure (use of single() cuts down on memory load)
PM = PAL_AMPM_setupPM('priorAlphaRange',single(priorAlphaRange),...
    'priorBetaRange',single(priorBetaRange),'priorGammaRange',single(gamma),...
    'priorLambdaRange',single(priorLambdaRange), 'numtrials',NumTrials, 'PF' , PF,...
    'stimRange',single(stimRange),'priorLambdaRange',single(priorLambdaRange),'marginalize',marginalize);

%trial loop
while PM.stop ~= 1

    drawnow
    clf

    response = rand(1) < PF(paramsGen, PM.xCurrent);    %simulate observer

%From here .......
    
    trackSuspend(length(PM.response)+1) = suspend;
    totals(find(single(stimRange)==PM.xCurrent)) = totals(find(single(stimRange)==PM.xCurrent)) + 1;

    if response == 1
        posterior = posterior.*vizLUT(:,:,:,find(PM.stimRange == PM.xCurrent));
    else
        posterior = posterior.*(1-vizLUT(:,:,:,find(PM.stimRange == PM.xCurrent)));
    end

    posterior = posterior./sum(sum(sum(posterior)));
    
    [maxim, I] = PAL_findMax(posterior);

    axes('units','pixels','position',[75 225 400 400]);
    
    %%%%%%%% (a): full 3-D posterior
    slice(a,b,l,PAL_Scale0to1(posterior)*64,[],[],lambdas);
    shading flat;
    if ~exist('OCTAVE_VERSION');
        alpha('color');
        am = linspace(.25,1,64);
        alphamap(am);
    end
    axis([min(alphas) max(alphas) min(betas) max(betas) min(lambdas) max(lambdas)])
    set(gca,'cameraposition',[-2.52 -1.57 0.1])
    set(gca,'xtick',[-.8:.4:.8]);
    set(gca,'ytick',log10([.1 .3 1 3 15]),'yticklabel',{'0.1','0.3','1','3','15'});
    set(gca,'xgrid', 'off','ygrid','off','zgrid','off')
    set(gca,'fontsize',8)
    xlabel('alpha','fontsize',12)
    ylabel('beta','fontsize',12)
    zlabel('lambda','fontsize',12)
    
    %%%%%%%% (b):
    
    if method == 1 || (method == 2 && sum(marginalize == 2))
        axes('units','pixels','position',[560 470 125 125]);
        toplot = squeeze(sum(sum(posterior(:,:,:),1),3))/max(squeeze(sum(sum(posterior(:,:,:),1),3)));
        plot(alphas,toplot,'k-')
        xlim = [min(alphas) max(alphas)];
        set(gca,'xlim',xlim)
        ylim = [0 1.2];
        set(gca,'ylim',ylim);
        set(gca,'xtick',[-.8:.4:.8],'ytick',[],'fontsize',8)
        xlabel('alpha','fontsize',10);
        text(xlim(1),ylim(2),'\beta, \lambda marginalized','fontsize',10,'verticalalignment','bottom')
        set(gca,'xdir','normal')
        set(gca,'ydir','normal')  
        
    else
        axes('units','pixels','position',[560 470 125 125]);
        xlim = [min(alphas) max(alphas)];
        set(gca,'xlim',xlim)
        ylim = [min(betas) max(betas)];
        set(gca,'ylim',ylim);
        set(gca,'ytick',log10([.1 .3 1 3 15]),'yticklabel',{'0.1','0.3','1','3','15'},'fontsize',8)
        set(gca,'xtick',[-.8:.4:.8],'fontsize',8)
        xlabel('alpha','fontsize',10);
        ylabel('beta','fontsize',10)
        set(gca,'tickdir','out')        
        if method == 0
            text(xlim(1),ylim(2),'\lambda fixed @ .03','fontsize',10,'verticalalignment','bottom')
            toplot = (squeeze(posterior(:,:,find(lambdas == 0.03)))./max(max(max(posterior))))*64;
        else
            text(xlim(1),ylim(2),'\lambda marginalized','fontsize',10,'verticalalignment','bottom')
            toplot = (squeeze(sum(posterior,3))./max(max(sum(posterior,3))))*64;
        end
        axes('units','pixels','position',[560 470 125 125]);        
        image(alphas, betas, toplot)
        set(gca,'xdir','normal')
        set(gca,'ydir','normal')
        set(gca,'ytick',[], 'xtick',[],'fontsize',6)
        
    end
    
    %%%%%% (c)
    if method == 0 || method == 1 || (method == 2 && sum(marginalize == 2))
        axes('units','pixels','position',[755 470 125 125]);
        if method == 0
            toplot = squeeze(sum(posterior(:,:,find(lambdas == 0.03)),1))/max(squeeze(sum(posterior(:,:,find(lambdas == 0.03)),1)));
            plot(alphas,toplot,'k-')
            xlim = [min(alphas) max(alphas)];
            set(gca,'xtick',[-.8:.4:.8],'ytick',[],'fontsize',8)
            xlabel('alpha','fontsize',10);
            text(xlim(1),ylim(2),'\lambda fixed, \beta marginalized','fontsize',10,'verticalalignment','bottom')
            set(gca,'xdir','normal')
            
        else
            toplot = squeeze(sum(sum(posterior(:,:,:),2),3))/max(squeeze(sum(sum(posterior(:,:,:),2),3)));
            plot(betas,toplot,'k-')
            xlim = [min(betas) max(betas)];
            set(gca,'xtick',log10([.1 .3 1 3 15]),'xticklabel',{'0.1','0.3','1','3','15'},'ytick',[],'fontsize',8)
            xlabel('beta','fontsize',10);
            text(xlim(2),ylim(2),'\alpha, \lambda marginalized','fontsize',10,'verticalalignment','bottom')
            set(gca,'xdir','reverse')
            
        end
        set(gca,'xlim',xlim)
        ylim = [0 1.2];
        set(gca,'ylim',ylim);
    else
        axes('units','pixels','position',[755 470 125 125]);
        xlim = [min(alphas) max(alphas)];
        set(gca,'xlim',xlim)
        ylim = [0 .1];
        set(gca,'ylim',ylim);
        set(gca,'ytick',[0:.02:.1],'fontsize',8)
        set(gca,'xtick',[-.8:.4:.8],'fontsize',8)
        xlabel('alpha','fontsize',10);
        ylabel('lambda','fontsize',10)
        set(gca,'tickdir','out')        
        text(xlim(1),ylim(2),'\beta marginalized','fontsize',10,'verticalalignment','bottom')
        toplot = (squeeze(sum(posterior,1))./max(max(sum(posterior,1))))'*64;
        axes('units','pixels','position',[755 470 125 125]);        
        image(alphas, lambdas, toplot)
        set(gca,'xdir','normal')
        set(gca,'ydir','normal')
        set(gca,'ytick',[], 'xtick',[],'fontsize',6)
        
    end

    %%%%%% (d)
    if method == 0 || method == 1 || (method == 2 && sum(marginalize == 2))
        axes('units','pixels','position',[950 470 125 125]);
        if method == 0
            toplot = squeeze(sum(posterior(:,:,find(lambdas == 0.03)),2))/max(squeeze(sum(posterior(:,:,find(lambdas == 0.03)),2)));
            plot(betas,toplot,'k-')
            xlim = [min(betas) max(betas)];
            set(gca,'xtick',log10([.1 .3 1 3 15]),'xticklabel',{'0.1','0.3','1','3','15'},'ytick',[],'fontsize',8)
            xlabel('beta','fontsize',10);
            text(xlim(2),ylim(2),'\lambda fixed, \alpha marginalized','fontsize',10,'verticalalignment','bottom')
            set(gca,'xdir','reverse')
        else
            toplot = squeeze(sum(sum(posterior(:,:,:),1),2))/max(squeeze(sum(sum(posterior(:,:,:),1),2)));
            plot(0:.01:.1,toplot,'k-')
            xlim = [0 .1];
            set(gca,'xtick',[0:.1:.1],'ytick',[],'fontsize',8)
            xlabel('lambda','fontsize',10);
            text(xlim(1),ylim(2),'\alpha, \beta marginalized','fontsize',10,'verticalalignment','bottom')
            set(gca,'xdir','normal')
        end
        set(gca,'xlim',xlim)
        ylim = [0 1.2];
        set(gca,'ylim',ylim);
        
    else
        axes('units','pixels','position',[950 470 125 125]);
        xlim = [min(betas) max(betas)];
        set(gca,'xlim',xlim)
        ylim = [0 .1];
        set(gca,'ylim',ylim);
        set(gca,'ytick',[0:.02:.1],'fontsize',8)
        set(gca,'xtick',log10([.1 .3 1 3 15]),'xticklabel',{'0.1','0.3','1','3','15'},'fontsize',8)
        xlabel('beta','fontsize',10);
        ylabel('lambda','fontsize',10)
        set(gca,'xdir','reverse')
        set(gca,'tickdir','out')        
        text(xlim(2),ylim(2),'\alpha marginalized','fontsize',10,'verticalalignment','bottom')
        toplot = (squeeze(sum(posterior,2))./max(max(sum(posterior,2))))'*64;
        axes('units','pixels','position',[950 470 125 125]);        
        image(betas, lambdas, toplot)
        set(gca,'xdir','reverse')
        set(gca,'ydir','normal')
        set(gca,'ytick',[], 'xtick',[],'fontsize',6)
    end

 %..... to here serves plotting only


 %Decide whether to suspend psi-marginal temporarily. Strategy effectively
 %draws WaitTime from exponential mass function resulting in constant
 %'hazard' (when in suspended mode, there is a constant probability (equal
 %to 1/WaitTime) of returning to psi-marginal mode on each trial).
 
    if PM.xCurrent == max(single(stimRange)) && AvoidConsecutive
        suspend = 1;
    end
    if suspend == 1
        suspend = rand(1) > 1./WaitTime;
    end            
    
 %update PM based on response
    PM = PAL_AMPM_updatePM(PM,response,'fixLapse',suspend);

    
 %Back to plotting.....

    %%(f)
    axes('units','pixels','position',[75 50 1000 125]);    
    t = 1:length(PM.x)-1;
    plot(t,PM.x(1:length(t)),'k');
    hold on
    plot(1:length(t),PM.threshold,'b-','LineWidth',2)
    plot(t(PM.response == 1 & trackSuspend == 0),PM.x(PM.response == 1 & trackSuspend == 0),'ko', 'MarkerFaceColor','k');
    plot(t(PM.response == 0 & trackSuspend == 0),PM.x(PM.response == 0 & trackSuspend == 0),'ko', 'MarkerFaceColor','w');
    plot(t(PM.response == 1 & trackSuspend == 1),PM.x(PM.response == 1 & trackSuspend == 1),'ro', 'MarkerFaceColor','r');
    plot(t(PM.response == 0 & trackSuspend == 1),PM.x(PM.response == 0 & trackSuspend == 1),'ro', 'MarkerFaceColor','w');
    for SR = 1:length(totals)
        if totals(SR) > 0
            plot(max(103,length(PM.x)+3),stimRange(SR),'ko','markerfacecolor','k','markersize',20*sqrt(totals(SR)./sum(totals)))
        end
    end
    minX = max(length(PM.x)-100, 0);
    maxX = max(length(PM.x)-100, 0)+105;
    minY = min(stimRange)-(max(stimRange)-min(stimRange))/10;
    maxY = max(stimRange)+(max(stimRange)-min(stimRange))/4;
    axis([minX maxX minY maxY])
    xlabel('Trial','fontsize',10);
    ylabel('log(intensity) (\itx\rm)','fontsize',10);
    plot(minX + 5, maxY - (maxY-minY)/12, 'ko', 'markerfacecolor','k')
    text(double(minX + 6), double(maxY - (maxY-minY)/12), 'correct')    
    plot(minX + 12, maxY - (maxY-minY)/12, 'ko', 'markerfacecolor','w')
    text(double(minX + 13), double(maxY - (maxY-minY)/12), 'incorrect') 
    if AvoidConsecutive
        plot(minX + 20, maxY - (maxY-minY)/12, 'ro', 'markerfacecolor','r')
        plot(minX + 22, maxY - (maxY-minY)/12, 'ro', 'markerfacecolor','w')
        text(double(minX + 23), double(maxY - (maxY-minY)/12), 'Selected while lapse rate was fixed at current ML estimate in attempt to avoid another ''free trial''') 
    end
       
    %%(e)
    axes('units','pixels','position',[560 240 515 175]);
    
    hold on
    [SL, NP, OON] = PAL_PFML_GroupTrialsbyX(PM.x(1:length(PM.x)-1),PM.response,ones(size(PM.response)));
    for SR = 1:length(SL(OON~=0))
        plot(SL(SR),NP(SR)/OON(SR),'ko','markerfacecolor','k','markersize',20*sqrt(OON(SR)./sum(OON)))
    end
    axis([min(stimRange)-(max(stimRange)-min(stimRange))/10 max(stimRange)+(max(stimRange)-min(stimRange))/10 0 1]);
    
    %crude(ish) ML fit
    [trash, MLindex] = PAL_findMax(PM.pdf);
    MLalpha = priorAlphaRange(MLindex(1));
    MLbeta = priorBetaRange(MLindex(2));
    if method == 0
        MLlambda = 0.03;
    else
        MLlambda = priorLambdaRange(MLindex(4));
    end
        
    plot([min(stimRange):.01:max(stimRange)],PF([0 1 gamma lambda],min(stimRange):.01:max(stimRange)),'k-','linewidth',2)
    plot([min(stimRange):.01:max(stimRange)],PF([MLalpha 10.^MLbeta gamma MLlambda],min(stimRange):.01:max(stimRange)),'r-','linewidth',2)
    plot([min(stimRange):.01:max(stimRange)],PF([PM.threshold(length(PM.threshold)) 10.^PM.slope(length(PM.threshold)) gamma PM.lapse(length(PM.threshold))],min(stimRange):.01:max(stimRange)),'b-','linewidth',2)
    
    xlabel('log(intensity) (\itx\rm)','fontsize',12);
    ylabel('\psi(\itx\rm; \alpha, \beta, \gamma, \lambda)','fontsize',12);
    text(min(stimRange)+(max(stimRange)-min(stimRange))/4, .25,'Bayes','color','b','Fontsize',14)
    text(min(stimRange)+(max(stimRange)-min(stimRange))/4, .1,'ML (crude)','color','r','Fontsize',14)
    set(gca,'xtick',[-1:.5:1]);
    set(gca,'ytick',[0:.25:1]);
    
    %inset
    axes('units','pixels','position',[955 255 100 100]);    
    hold on
    plot([min(stimRange):.01:max(stimRange)],PF([0 1],min(stimRange):.01:max(stimRange)),'k-','linewidth',2)
    plot([min(stimRange):.01:max(stimRange)],PF([MLalpha 10.^MLbeta],min(stimRange):.01:max(stimRange)),'r-','linewidth',2)
    plot([min(stimRange):.01:max(stimRange)],PF([PM.threshold(length(PM.threshold)) 10.^PM.slope(length(PM.threshold))],min(stimRange):.01:max(stimRange)),'b-','linewidth',2)
    set(gca,'xtick',[],'ytick',[0 1],'Fontsize',10);
    ylabel('F(\itx\rm)')
    
    drawnow
    
    %Allow user input to step through (space bar) or automate (g)
    if ~exist('OCTAVE_VERSION');
        cc = get(gcf,'currentcharacter');
        if cc == ' ';
            pause
        end   
        if length(PM.response) == 1
            pause
        end
    end
end