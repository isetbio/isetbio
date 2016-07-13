%
%PAL_MLDS_Demo  Demonstrates (1) generating stimulus lists for a Maximum
%Likelihood Difference Scaling experiment involving stimulus pairs, triads 
%or quads (2) fitting the results of such an experiment using a Maximum
%Likelihood criterion and (3) estimating standard errors on parameters
%using bootstrap simulations.
%
%Demonstrates usage of Palamedes functions:
%-PAL_MLDS_GenerateStimList
%-PAL_MLDS_SimulateObserver
%-PAL_MLDS_GroupTrialsbyX
%-PAL_MLDS_Fit
%-PAL_MLDS_Bootstrap
%secondary:
%-PAL_Logistic
%-PAL_Scale0to1
%
%More information on any of these functions may be found by typing
%help followed by the name of the function. e.g., help PAL_MLDS_Fit
%
%FK & NP (September 2009)

clear all;

ptq = input('Type 2 for pairs, 3 for triads, 4 for quadruples: ');
NumLevels = input('Number of stimulus levels (e.g., 10): ');
NumRepeats = input('Number of repetitions per comparison (e.g., 40): ');
SDnoise = input('Value for the noise standard deviation (e.g., .3): ');

if exist('OCTAVE_VERSION');
    fprintf('\nUnder Octave, Figure may not render exactly as intended. Visit\n');
    fprintf('www.palamedestoolbox.org/demosfiguregallery.html to see figure\n');
    fprintf('as intended.\n\n');
end

disp('Working. Please wait.');

%set options for Nelder-Mead simplex search (optional)
options = PAL_minimize('options');
options.MaxFunEvals = 2000;         %allow more function evaluations
options.MaxIter = 2000;             %allow more iterations
options.TolFun = 1e-09;             %require higher precision on LL
options.TolX = 1e-09;               %require higher precision on parameter
                                    %estimates

%Simulated observer's characteristics
PF = @PAL_Logistic;         %Shape
alpha=(NumLevels+1)/2;      %Threshold
beta=10/NumLevels;          %Slope
params = [alpha beta 0 0];  %Threshold Slope Guess rate=0 Lapse rate=0
PsiValuesGen = PF(params,[1:NumLevels]); %psi values for simulated observer

MaxDiff = 3;  %Stimuli (pairs: [i j], triads: [i j k], quads: [i j k l],  
              %i < j < k < l) will be generated for which:
              %pairs: j - i will be at most MaxDiff
              %triads: abs([k - j]-[j - i]) will be at most MaxDiff
              %quads: abs([l - k]-[j - i]) will be at most MaxDiff

B = 400;      %Number of bootstrap simulations to perform for SEs

%Generate stimulus list
StimList = PAL_MLDS_GenerateStimList(ptq, NumLevels, MaxDiff, NumRepeats);

%Each entry of StimList corresponds to single trial
OutOfNum = ones(1,size(StimList,1));

%Simulate observer
response = PAL_MLDS_SimulateObserver(StimList, OutOfNum, PsiValuesGen, ...
    SDnoise);

%The following groups identical trials together
%Before this line 'StimList' will have as many rows as there are trials, 
%after this line, 'StimList' will have as many rows as there are unique 
%stimuli. NumPos will contain the number of positive responses at each 
%trial type, OutOfNum will contain the total number of trials at which the 
%trial type was presented
[StimList NumPos OutOfNum] = PAL_MLDS_GroupTrialsbyX(StimList, response,...
    OutOfNum);    

%Initial search values
PsiValuesStart=PAL_Scale0to1(PsiValuesGen);

%Perform fit
[PsiValues SDnoise LL exitflag output] = PAL_MLDS_Fit(StimList, NumPos, ...
    OutOfNum, PsiValuesStart, SDnoise, 'SearchOptions', options);

%Or just:
%[PsiValues SDnoise LL exitflag output] = PAL_MLDS_Fit(StimList, NumPos, ...
%    OutOfNum, PsiValuesStart, SDnoise);

%Perform bootstrap to find SEs
maxTries = 4;   %try fits up to 4 times
rangeTries = .1;  %range of random jitter on initial search values on 
                  %retried fits. or: rangeTries = .1.*ones(1,NumLevels-1);

[SE_PsiValues SE_SDnoise] = PAL_MLDS_Bootstrap(StimList,OutOfNum,...
    PsiValues,SDnoise,B, 'SearchOptions', options, 'maxTries', maxTries,...
    'RangeTries',rangeTries);

%calculate total number of trials
NumTrials=size(StimList,1).*NumRepeats;
message = ['Number of trials: ', num2str(NumTrials)];
disp(message);

%plot results
figure('name','Maximum Likelihood Difference Scaling');
set(gca, 'fontsize',18);
plot(1:NumLevels, PsiValues, 'k-s','Markersize',8);
for i = 2:length(SE_PsiValues)-1
    line([i i],[PsiValues(i)-SE_PsiValues(i) PsiValues(i)+...
        SE_PsiValues(i)], 'color','k');
end
axis ([0 NumLevels+1 -0.2 1.2]);
hold on

%Add plot of function used to generate input responses
StimLevelsGenPlot=[1:0.1:NumLevels];
PsiValuesGenPlot = PF(params,StimLevelsGenPlot);  
plot(StimLevelsGenPlot, PAL_Scale0to1(PsiValuesGenPlot), '-','LineWidth',2,'color',[0 .7 0]);
if exist('OCTAVE_VERSION')
    plot(1,1.12,'ks','MarkerSize',8);
    line([.6 1.4],[1.12 1.12],'color','k');
    line([.6 1.4],[1.03 1.03],'color',[0 .7 0],'linewidth',2);
    text(1.5,1.12,'Parameter estimates','fontsize',12);
    text(1.5,1.03,'Response generating function','fontsize',12);
    rectangle('Position',[.3 .98 6.8 .19])
else
    [h o] = legend('Parameter estimates','Response generating function',2);
    set(h,'Interpreter','none','fontsize',12,'Location','Northwest');
    set(o(5),'color',[0 .7 0],'LineWidth',2);
end
% Add some labeling
xlabel('Stimulus magnitude','FontSize',18);
ylabel('Perceptual magnitude','FontSize',18);
