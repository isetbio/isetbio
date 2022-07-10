function testSubjectID = subjectWithRankInEye(obj, subjectRankOrder, whichEye)

    % Ensure we have a valid eye specification
    assert(ismember(whichEye, RetinaToVisualFieldTransformer.validEyes), ...
        'Wrong eye specification: ''%s''.', whichEye);

    switch (obj.ZernikeDataBase)
        case RetinaToVisualFieldTransformer.Artal
            rankedSujectIDs = ArtalOptics.constants.subjectRanking(whichEye);
            testSubjectID = rankedSujectIDs(subjectRankOrder);
        case RetinaToVisualFieldTransformer.Polans
            if (~strcmp(whichEye, RetinaToVisualFieldTransformer.rightEye))
                error('Polans measurements exist only for the right eye.');
            end
            rankedSujectIDs = PolansOptics.constants.subjectRanking();
            testSubjectID = rankedSujectIDs(subjectRankOrder);

        otherwise
            error('Unknown zernike database: ''%ss'.', obj.ZernikeDataBase);
    end

end