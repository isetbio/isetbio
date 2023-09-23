function [stimulusChromaticity, responsesFileName] = ...
    chooseStimulusChromaticityForMosaicResponsesAndUpdateFileName(...
    responsesFileName, identifierString)

    stimulusChromaticity = 'invalid';
    validChromaticityChoices = {'achromatic', 'Lcone isolating', 'Mcone isolating', 'Scone isolating'};

    while (~ismember(stimulusChromaticity, validChromaticityChoices)) 
           fprintf('\nChoose stimulus chromaticity \n');
           fprintf('[A]  - Achromatic (default) \n');
           fprintf('[L]  - L-cone isolating \n');
           fprintf('[M]  - M-cone isolating \n');
           fprintf('[S]  - S-cone isolating \n');
           chromaChoice = input('Stimulus chromaticity: ', 's');

           stimulusChromaticity = 'invalid';
           switch (upper(chromaChoice))
                case 'A'
                    stimulusChromaticity = 'achromatic';
                    % We do not add a chromaticity in the filename if achromatic stimulus was employed - this is the default
                case 'L'
                    stimulusChromaticity = 'Lcone isolating';
                    responsesFileName = strrep(responsesFileName, identifierString, sprintf('LconeIsolating%s', identifierString));
                case 'M'
                    stimulusChromaticity = 'Mcone isolating';
                    responsesFileName = strrep(responsesFileName, identifierString, sprintf('MconeIsolating%s', identifierString));
                case 'S'
                    stimulusChromaticity = 'Scone isolating';
                    responsesFileName = strrep(responsesFileName, identifierString, sprintf('SconeIsolating%s', identifierString));
           end % switch

    end % while
end

