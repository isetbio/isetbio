%
%PAL_AMUD_Demo  Demonstrates use of Palamedes routines to implement an
%'up/down' adaptive procedure.
%
%Demonstrates usage of Palamedes functions:
%-PAL_AMUD_setupUD
%-PAL_AMUD_updateUD
%-PAL_AMUD_analyzeUD
%secondary:
%PAL_Gumbel
%PAL_PFML_Fit
%
%More information on any of these functions may be found by typing
%help followed by the name of the function. e.g., help PAL_AMUD_setupUD
%
%NP (September 2009)

clear all;  %Clear all existing variables from memory

%Simulated observer's characteristics
PF = @PAL_Gumbel;
trueParams = [0 20 .5 0.01];

%Set up up/down procedure:
up = 1;                     %increase after 1 wrong
down = 3;                   %decrease after 3 consecutive right
StepSizeDown = 0.05;        
StepSizeUp = 0.05;
stopcriterion = 'trials';   
stoprule = 50;
startvalue = 0.3;           %intensity on first trial

UD = PAL_AMUD_setupUD('up',up,'down',down);
UD = PAL_AMUD_setupUD(UD,'StepSizeDown',StepSizeDown,'StepSizeUp', ...
    StepSizeUp,'stopcriterion',stopcriterion,'stoprule',stoprule, ...
    'startvalue',startvalue);

%Determine and display targetd proportion correct and stimulus intensity
targetP = (StepSizeUp./(StepSizeUp+StepSizeDown)).^(1./down);
message = sprintf('\rTargeted proportion correct: %6.4f',targetP);
disp(message);
targetX = PAL_Gumbel(trueParams, targetP,'inverse');
message = sprintf('Targeted stimulus intensity given simulated observer');
message = strcat(message,sprintf(': %6.4f',targetX));
disp(message);

%Trial loop

while ~UD.stop

    %Present trial here at stimulus intensity UD.xCurrent and collect
    %response
    %Here we simulate a response instead (0: incorrect, 1: correct)
    response = rand(1) < PF(trueParams, UD.xCurrent);   
    
    UD = PAL_AMUD_updateUD(UD, response); %update UD structure
    
end

%Threshold estimate as mean of all but the first three reversal points
Mean = PAL_AMUD_analyzeUD(UD, 'reversals', max(UD.reversal)-3);
message = sprintf('\rThreshold estimate as mean of all but last three');
message = strcat(message,sprintf(' reversals: %6.4f', Mean));
disp(message);

%Threshold estimate found by fitting Gumbel
params = PAL_PFML_Fit(UD.x, UD.response, ones(1,length(UD.x)), ...
    [0 50 .5 .01], [1 0 0 0], @PAL_Gumbel);
message = sprintf('Threshold estimate as alpha of fitted Gumbel: %6.4f'...
    , params(1));
disp(message);

%Create simple plot:
t = 1:length(UD.x);
figure('name','Up/Down Adaptive Procedure');
plot(t,UD.x,'k');
hold on;
plot(t(UD.response == 1),UD.x(UD.response == 1),'ko', 'MarkerFaceColor','k');
plot(t(UD.response == 0),UD.x(UD.response == 0),'ko', 'MarkerFaceColor','w');
set(gca,'FontSize',16);
axis([0 max(t)+1 min(UD.x)-(max(UD.x)-min(UD.x))/10 max(UD.x)+(max(UD.x)-min(UD.x))/10]);
line([1 length(UD.x)], [targetX targetX],'linewidth', 2, 'linestyle', '--', 'color','k');
xlabel('Trial');
ylabel('Stimulus Intensity');