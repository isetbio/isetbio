function [stimulusChromaticity, coneMosaicSTFresponsesFileName] = ...
    chooseStimulusChromaticityForInputConeMosaicSTFresponses(coneMosaicSTFresponsesFileName)

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
                case 'L'
                    stimulusChromaticity = 'Lcone isolating';
                    coneMosaicSTFresponsesFileName = strrep(coneMosaicSTFresponsesFileName, 'coneMosaicSTFresponses', 'coneMosaicLconeIsolatingSTFresponses');
                case 'M'
                    stimulusChromaticity = 'Mcone isolating';
                    coneMosaicSTFresponsesFileName = strrep(coneMosaicSTFresponsesFileName, 'coneMosaicSTFresponses', 'coneMosaicMconeIsolatingSTFresponses');
                case 'S'
                    stimulusChromaticity = 'Scone isolating';
                    coneMosaicSTFresponsesFileName = strrep(coneMosaicSTFresponsesFileName, 'coneMosaicSTFresponses', 'coneMosaicSconeIsolatingSTFresponses');
           end % switch

    end % while
end

