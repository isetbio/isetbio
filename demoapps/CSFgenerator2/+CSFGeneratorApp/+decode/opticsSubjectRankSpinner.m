function opticsSubjectRankSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.opticsSubjectRankSpinner.Value = value;
        case 'sliderToValue'
            switch (app.opticsParams.subjectDataset)
                case 'Polans2015'
                    maxSubjectRank = app.opticsParams.maxPolansSubjectRank;
                case 'Artal2012'
                    maxSubjectRank = app.opticsParams.maxArtalSubjectRank;
                otherwise
                    error('Unknown Zernike data set: ''%s''.', app.opticsParams.subjectDataset)
            end
            % Change the range of available subject ranks
            app.opticsSubjectRankSpinner.Limits = [1 maxSubjectRank];
            % Change the value of the subject rank
            app.opticsParams.subjectRank = min([value maxSubjectRank]);
            app.opticsSubjectRankSpinner.Value = app.opticsParams.subjectRank;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end
