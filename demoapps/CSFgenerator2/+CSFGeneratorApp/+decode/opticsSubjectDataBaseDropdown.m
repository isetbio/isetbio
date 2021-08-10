function opticsSubjectDataBaseDropdown(app, direction, value)
    
    switch direction
        case 'valueToSlider'
            app.opticsSubjectDataBaseDropDown.Value = value;
        case 'sliderToValue'
            app.opticsParams.subjectDataset = value;
            switch (value)
                case 'Polans2015'
                    maxSubjectRank = app.opticsParams.maxPolansSubjectRank;
                    
                    % Only right eye data exist, so disable eye switch
                    app.roiEyeSwitch.Enable = 'off';
                    
                    % We were in left eye, so we need to switch to right
                    % eye. This means new processing pipeline
                    if strcmp(app.coneMosaicParams.whichEye, 'left eye')  
                        app.roiParams.whichEye = 'right eye';
                        app.roiEyeSwitch.Value = 'right eye';
                        app.coneMosaicParams.whichEye = 'right eye';
                        
                        CSFGeneratorApp.generate.processingPipeline(app, 'Generating new processing pipeline');
                        CSFGeneratorApp.render.coneMosaicView(app, 'update');
                        CSFGeneratorApp.render.roiView(app, 'update');
                    end
                    
                case 'Artal2012'
                    maxSubjectRank = app.opticsParams.maxArtalSubjectRank;
                    % Data for both eyes exist, so enable eye switch
                    app.roiEyeSwitch.Enable = 'on';
                otherwise
                    error('Unknown Zernike data set: ''%s''.', app.opticsParams.subjectDataset)
            end
            % Change the range of available subject ranks
            previousRank = app.opticsParams.subjectRank;
            app.opticsSubjectRankSpinner.Value = 1;
            app.opticsSubjectRankSpinner.Limits = [1 maxSubjectRank];
            % Change the value of the subject rank
            app.opticsSubjectRankSpinner.Value = min([previousRank maxSubjectRank]);
            app.opticsParams.subjectRank = app.opticsSubjectRankSpinner.Value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
    
    
end
