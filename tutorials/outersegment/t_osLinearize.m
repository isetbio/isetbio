% t_osLinearize
%
% Illustrates thefoveal and peripheral cone outer segment tempral impulse
% responses using the osBioPhys object at different mean intensities.
%
% A single cone is created, and its absorption time course is set to have an
% impulse at the specified time step. Biophysical outer segments are created,
% one with foveal dynamics and one with peripheral dynamics.
%
% The goal is to create a lookup table of the impulse response at different
% intensities.
%
% 11/2016 JRG (c) isetbio team

%% 
ieInit;

%% Peripheral cone dynamics

% Mean background isomerization in R*/sec
% We will use a range of these, but here are some typical values
% Changing this value only changes the starting current, which we subtract
% out, so this doesn't affect the impulse response.
backgroundMean = 1000;   %[10 100 500 1000 2000 5000 10000];

% Mean of isomerization stimulus in R*/sec (scaled by time step to be
% placed in bins of photons). The mean rate affects the magnitude of
% the impulse response due to adaptation.
meanIsoArray = [100 500 1000 2000 5000 10000]*timeStep;

% Set up parameters for stimulus.
nSamples = 8000;        % Number of time samples 
timeStep = 1e-4;        % time step

% Extra photons at the flash time bin
% Units are R*
flashIntens   = 10;      %
flashTimeStep = 5000;    % Warm up period, then flash

vcNewGraphWin; hold on

for meanInd = 1:length(meanIsoArray)
    for bgInd = 1:length(backgroundMean)
        
        % Build the stimulus
        stimulus = ones(1,nSamples)*backgroundMean(bgInd);
        
        % Set mean background for start of current output
        bgIsomerizations = backgroundMean(bgInd);
        
        % flash intensity in R*/cone/sec (maintained for 1 bin only)
        meanIntens  = meanIsoArray(meanInd);      
        
        % Create stimulus.
        stimulus = meanIntens*ones(nSamples, 1);
        stimulusFlat = reshape(stimulus, [1 1 nSamples]);
        
        % stimulus(1) = flashIntens*timeStep;
        stimulus(flashTimeStep) = flashIntens;
        stimulus = reshape(stimulus, [1 1 nSamples]);
        
        % Generate the cone mosaic for the impulse stimulus
        osCM = osBioPhys('osType',false);            % peripheral (fast) cone dynamics
        osCM.set('noise flag',0);
        cm = coneMosaic('os',osCM,'pattern', 2); % a single cone
        cm.integrationTime = timeStep;
        cm.os.timeStep = timeStep;
        
        cm.absorptions  = stimulus;
        
        % Compute outer segment currents.
        cm.computeCurrent('bgR',bgIsomerizations);
        currentScaled = (cm.current);% - cm.current(flashTimeStep);
        
        % Generate the cone mosaic for the flat stimulus
        osCMflat = osBioPhys('osType',false);            % peripheral (fast) cone dynamics
        osCMflat.set('noise flag',0);
        cmFlat = coneMosaic('os',osCMflat ,'pattern', 2); % a single cone
        cmFlat.integrationTime = timeStep;
        cmFlat.os.timeStep = timeStep;
        
        cmFlat.absorptions  = stimulusFlat;
        
        % Compute outer segment currents.
        cmFlat.computeCurrent('bgR',bgIsomerizations);
        currentScaledFlat = (cmFlat.current);% - cm.current(flashTimeStep);
        
        % Plot the impulse responses for comparison.
        
        subplot(1,2,1); hold on;
        tme = (1:nSamples)*timeStep;
        plot(tme,squeeze(currentScaled),'LineWidth',3);
        grid on; % legend('Peripheral','Foveal');
        % axis([0 0.2 min((current2Scaled./max(current2Scaled(:)))) 1]);
        xlabel('Time (sec)','FontSize',14);
        ylabel('pA','FontSize',14);
        set(gca,'fontsize',14);
        title(sprintf('Peripheral (fast)\nImpulse response as a function of mean R*/sec'));
        
        subplot(1,2,2); hold on;
        tme = (1:nSamples)*timeStep;
        plot(tme-flashTimeStep*timeStep,squeeze(currentScaled)- squeeze(currentScaledFlat),'LineWidth',3);
        grid on; % legend('Peripheral','Foveal');
        % axis([0 0.2 min((current2Scaled./max(current2Scaled(:)))) 1]);
        axis([0 nSamples*timeStep-flashTimeStep*timeStep -.2 2]);
        xlabel('Time (sec)','FontSize',14);
        ylabel('pA','FontSize',14);
        set(gca,'fontsize',14);
        % title('Impulse response in the dark
        
        title(sprintf('Peripheral (fast)\nIR as a function of mean R*/sec, rezeroed'));
        
    end
end
% legend(num2str(meanIsoArray));
legend('100', '500', '1000', '2000', '5000', '10000')
