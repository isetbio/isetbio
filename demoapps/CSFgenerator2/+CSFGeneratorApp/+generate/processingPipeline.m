function processingPipeline(app, headerMessage)
    
     % Determine whether we will show any more progress bars
     %app.showsProgressBar = getpref('ISET','waitbar');
            
     % Start the parpool
     dialogBox = uiprogressdlg(app.mainView,'Title',headerMessage,...
                'Message',sprintf('This app requires MATLAB 2020b or later, and that (i) isetbio, (ii) ISETBioCSFGenerator, and (iii) mQuestPlus, available at: \n\t-https://github.com/isetbio/isetbio\n\t-https://github.com/isetbio/ISETBioCSFGenerator.git\n\t-https://github.com/BrainardLab/mQUESTPlus.git \nare on the user''s path.\n\nWaiting for parallel pool engine to wake up...'));
     
     dialogBox.Value = 0.2; 
     
     poolobj = gcp('nocreate'); 
     if isempty(poolobj)
         pause(.1)
         parpool('local');
     end
                
    dialogBox.Value = 0.5; 
    dialogBox.Message = 'Generating cone mosaic. Please wait ...';
    CSFGeneratorApp.generate.coneMosaic(app, dialogBox);
    
    
    dialogBox.Value = 0.8; 
    dialogBox.Message = 'Generating optics. Please wait ...';        
    CSFGeneratorApp.generate.optics(app, dialogBox);
    
    dialogBox.Value = 0.9; 
    dialogBox.Message = 'Generating stimulus. Please wait ...';
    CSFGeneratorApp.generate.stimulus(app, dialogBox);
    
    close(dialogBox);
end

