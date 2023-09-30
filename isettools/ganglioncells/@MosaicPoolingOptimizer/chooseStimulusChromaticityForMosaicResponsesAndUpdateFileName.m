function [stimulusChromaticity, responsesFileName] = ...
    chooseStimulusChromaticityForMosaicResponsesAndUpdateFileName(...
    responsesFileName, identifierString, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('doNotQueryUserInsteadEmployThisStimulusChromaticity', '', @ischar);
    p.parse(varargin{:});

    doNotQueryUserInsteadEmployThisStimulusChromaticity = p.Results.doNotQueryUserInsteadEmployThisStimulusChromaticity;

    if (~isempty(doNotQueryUserInsteadEmployThisStimulusChromaticity))
        stimulusChromaticity = doNotQueryUserInsteadEmployThisStimulusChromaticity;
        responsesFileName = updateResponsesFileName(stimulusChromaticity, responsesFileName, identifierString);
        return;
    end

    stimulusChromaticity = 'invalid';
    
    while (~ismember(stimulusChromaticity, MosaicPoolingOptimizer.validChromaticityChoices)) 
           fprintf('\nChoose stimulus chromaticity \n');
           fprintf('[A]  - Achromatic (default) \n');
           fprintf('[L]  - L-cone isolating \n');
           fprintf('[M]  - M-cone isolating \n');
           fprintf('[S]  - S-cone isolating \n');
           chromaChoice = input('Stimulus chromaticity: ', 's');

           switch (upper(chromaChoice))
                case 'A'
                    stimulusChromaticity = 'achromatic';
                    % We do not add a chromaticity in the filename if achromatic stimulus was employed - this is the default
                case 'L'
                    stimulusChromaticity = 'Lcone isolating';
                    responsesFileName = updateResponsesFileName(stimulusChromaticity, responsesFileName, identifierString);
                case 'M'
                    stimulusChromaticity = 'Mcone isolating';
                    responsesFileName = updateResponsesFileName(stimulusChromaticity, responsesFileName, identifierString);
                case 'S'
                    stimulusChromaticity = 'Scone isolating';
                    responsesFileName = updateResponsesFileName(stimulusChromaticity, responsesFileName, identifierString);
           end % switch
    end % while
end


function responsesFileName = updateResponsesFileName(stimulusChromaticity, responsesFileName, identifierString)
    
    switch (stimulusChromaticity)
        case 'achromatic'
            % We do not add a chromaticity in the filename if achromatic stimulus was employed - this is the default
        case 'Lcone isolating'
            responsesFileName = strrep(responsesFileName, identifierString, sprintf('LconeIsolating%s', identifierString));
        case 'Mcone isolating'
            responsesFileName = strrep(responsesFileName, identifierString, sprintf('MconeIsolating%s', identifierString));
        case 'Scone isolating'
            responsesFileName = strrep(responsesFileName, identifierString, sprintf('SconeIsolating%s', identifierString));
        otherwise
            error('Unknown stimulus chromaticity: ''%s''.', stimulusChromaticity);

     end % switch

end
