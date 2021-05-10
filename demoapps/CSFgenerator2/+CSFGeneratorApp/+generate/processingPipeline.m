function processingPipeline(app)
    
    dialogBox = uiprogressdlg(app.mainView,'Title','Please Wait',...
                    'Message','');
                
    dialogBox.Value = 0.3; 
    dialogBox.Message = 'Generating cone mosaic. Please wait ...';
    CSFGeneratorApp.generate.coneMosaic(app);
    
    
    dialogBox.Value = 0.6; 
    dialogBox.Message = 'Generating optics. Please wait ...';        
    CSFGeneratorApp.generate.optics(app);
    
    close(dialogBox);
end

