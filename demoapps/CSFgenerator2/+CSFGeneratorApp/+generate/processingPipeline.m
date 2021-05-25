function processingPipeline(app, headerMessage, varargin)
    
     p = inputParser;
     p.addParameter('doNotRegenerateConeMosaic', false, @islogical);
     p.addParameter('doNotRegenerateOptics', false, @islogical);
     p.addParameter('doNotRegenerateStimulus', false, @islogical);
     p.parse(varargin{:});
     
     regenerateConeMosaic = ~p.Results.doNotRegenerateConeMosaic;
     regenerateOptics = ~p.Results.doNotRegenerateOptics;
     regenerateStimulus = ~p.Results.doNotRegenerateStimulus;
     
     % Determine whether we will show any more progress bars
     %app.showsProgressBar = getpref('ISET','waitbar');
            
     % Start the parpool
     dialogBox = uiprogressdlg(app.mainView,'Title',headerMessage,...
                'Message',sprintf('This app requires MATLAB 2020b or later, and that (i) isetbio, (ii) ISETBioCSFGenerator, and (iii) mQuestPlus, available at: \n\t-https://github.com/isetbio/isetbio\n\t-https://github.com/isetbio/ISETBioCSFGenerator.git\n\t-https://github.com/BrainardLab/mQUESTPlus.git \nare on the user''s path.\n\nWaiting for parallel pool engine to wake up...'));
     
     dialogBox.Value = 0.2; 
     
     poolobj = gcp('nocreate'); 
     if isempty(poolobj)
         drawnow;
         pause(0.5)
         parpool('local');
     end
         
    if (regenerateConeMosaic)
        % Generate the cone mosaic
        dialogBox.Value = 0.5; 
        dialogBox.Message = 'Generating cone mosaic. Please wait ...';
        CSFGeneratorApp.generate.coneMosaic(app, dialogBox);
    end
    
    if (regenerateOptics)
        % Generate the optics
        dialogBox.Value = 0.6; 
        dialogBox.Message = 'Generating optics. Please wait ...';        
        CSFGeneratorApp.generate.optics(app, dialogBox);
    end
    
    if (regenerateStimulus)
        % Compute the stimulus
        dialogBox.Value = 0.7; 
        dialogBox.Message = 'Generating stimulus. Please wait ...';
        CSFGeneratorApp.generate.stimulus(app, dialogBox);
    end
    
    % Compute the cone mosaic activation
    dialogBox.Value = 0.8; 
    dialogBox.Message = 'Computing cone mosaic response. Please wait ...';
    CSFGeneratorApp.compute.coneMosaicActivation(app, dialogBox);
    
    close(dialogBox);
end

