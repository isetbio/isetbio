%
%PAL_AMRF_Demo  Demonstrates use of Palamedes routines to implement a
%'running fit' adaptive procedure (e.g., best PEST, QUEST).
%
%Halfway (after 25 trials) run is interrupted, intermediate results
%are stored to disk, RAM memory is cleared, intermediate results are
%reloaded from disk, and run continues (as if observer runs 2 sessions on
%separate occassions). The posterior from first session serves as prior to
%second session.
%
%Demonstrates usage of Palamedes functions:
%-PAL_AMRF_setupRF
%-PAL_AMRF_updateRF
%secondary:
%PAL_Gumbel
%PAL_pdfNormal
%
%More information on any of these functions may be found by typing
%help followed by the name of the function. e.g., help PAL_AMRF_setupRF
%
%NP (September 2009)

clear all;  %Clear all existing variables from memory

%Simulated observer's characteristics
PFsimul = @PAL_Gumbel;
trueParams = [0 2 0.5 0.01];

%Set up running fit procedure:

%Define prior
alphas = [-2:.01:2];
prior = PAL_pdfNormal(alphas,0,2); %Gaussian

%Termination rule
stopcriterion = 'trials';
stoprule = 25;

%Function to be fitted during procedure
PFfit = @PAL_Gumbel;    %Shape to be assumed
beta = 2;               %Slope to be assumed
lambda  = 0.01;         %Lapse rate to be assumed
meanmode = 'mean';      %Use mean of posterior as placement rule

%set up procedure
RF = PAL_AMRF_setupRF('priorAlphaRange', alphas, 'prior', prior,...
    'stopcriterion',stopcriterion,'stoprule',stoprule,'beta',beta,...
    'lambda',lambda,'PF',PFfit,'meanmode',meanmode);

%Trial loop

while RF.stop ~= 1
    
    %Present trial here at stimulus intensity UD.xCurrent and collect
    %response
    %Here we simulate a response instead (0: incorrect, 1: correct)    
    amplitude = RF.xCurrent;
    response = rand(1) < PFsimul(trueParams,amplitude);    
    RF = PAL_AMRF_updateRF(RF, amplitude, response);
    
end

save savedRF RF  %by saving RF structure all data are saved
clear all   %this clears all variables from RAM memory (as if you were to  
            %shut down Matlab and computer).
load savedRF     %load RF structure to pick up where you left off

%Simulated observer's characteristics (these were deleted by 'clear all')
PFsimul = @PAL_Gumbel;
trueParams = [0 2 0.5 0.01];

%Change stop rule to 50 total trials and set RF.stop back to 0.
%In effect, this uses posterior of previous session as prior for new
%session.
RF = PAL_AMRF_setupRF(RF,'stoprule',50);
RF.stop = 0;

%Trial loop

while RF.stop ~= 1
    
    amplitude = RF.xCurrent;
    response = rand(1) < PFsimul(trueParams,amplitude);    
    RF = PAL_AMRF_updateRF(RF, amplitude, response);
    
end

%Print summary of results to screen
message = sprintf('\rThreshold estimate as mode of posterior: %6.4f'...
    ,RF.mode);
disp(message);
message = sprintf('Threshold estimate as mean of posterior: %6.4f'...
    ,RF.mean);
disp(message);
message = sprintf('Threshold standard error as sd of posterior: %6.4f'...
    ,RF.sd);
disp(message);

message = sprintf('\rThreshold estimate as mode of posterior using unifo');
message = strcat(message,sprintf('rm prior (i.e., ML estimate): %6.4f' ,...
    RF.modeUniformPrior));
disp(message);
message = sprintf('Threshold estimate as mean of posterior using unifo');
message = strcat(message,sprintf('rm prior: %6.4f' ,...
    RF.meanUniformPrior));
disp(message);
message = sprintf('Threshold standard error as sd of posterior using uni');
message = strcat(message,sprintf('form prior: %6.4f' ,...
    RF.sdUniformPrior));
disp(message);

%Create simple plot:
t = 1:length(RF.x);
figure('name','Running Fit Adaptive Procedure');
plot(t,RF.x,'k');
hold on;
plot(t(RF.response == 1),RF.x(RF.response == 1),'ko', ...
    'MarkerFaceColor','k');
plot(t(RF.response == 0),RF.x(RF.response == 0),'ko', ...
    'MarkerFaceColor','w');
set(gca,'FontSize',16);
axis([0 max(t)+1 min(RF.x)-(max(RF.x)-min(RF.x))/10 ...
    max(RF.x)+(max(RF.x)-min(RF.x))/10]);
line([1 length(RF.x)], [trueParams(1) trueParams(1)],'linewidth', 2, ...
    'linestyle', '--', 'color','k');
xlabel('Trial');
ylabel('Stimulus Intensity');