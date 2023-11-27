function [stimulusChromaticity, coneFundamentalsOptimizedForStimPosition, responsesFileName] = ...
    chooseStimulusChromaticityForMosaicResponsesAndUpdateFileName(...
    responsesFileName, identifierString, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('doNotQueryUserButUseThisStimulusChromaticityAndConeFundamentals', [], @(x)(isempty(x)||(isstruct(x))));
    p.parse(varargin{:});

    doNotQueryUserButUseThisStimulusChromaticityAndConeFundamentals = p.Results.doNotQueryUserButUseThisStimulusChromaticityAndConeFundamentals;

    if (~isempty(doNotQueryUserButUseThisStimulusChromaticityAndConeFundamentals))
        stimulusChromaticity = doNotQueryUserButUseThisStimulusChromaticityAndConeFundamentals.stimulusChromaticity;
        if (~isempty(doNotQueryUserButUseThisStimulusChromaticityAndConeFundamentals.coneFundamentalsOptimizedForStimPosition))
            coneFundamentalsOptimizedForStimPosition = doNotQueryUserButUseThisStimulusChromaticityAndConeFundamentals.coneFundamentalsOptimizedForStimPosition;
        else
            coneFundamentalsOptimizedForStimPosition = queryUserWhetherToEmployOptimizedConeFundamentals();
        end
        
        responsesFileName = updateResponsesFileName(stimulusChromaticity, coneFundamentalsOptimizedForStimPosition, responsesFileName, identifierString);
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
                    coneFundamentalsOptimizedForStimPosition = queryUserWhetherToEmployOptimizedConeFundamentals();
                    % We do not add a chromaticity in the filename if achromatic stimulus was employed - this is the default
                case 'L'
                    coneFundamentalsOptimizedForStimPosition = queryUserWhetherToEmployOptimizedConeFundamentals();
                    stimulusChromaticity = 'Lcone isolating';
                    responsesFileName = updateResponsesFileName(stimulusChromaticity, coneFundamentalsOptimizedForStimPosition, responsesFileName, identifierString);
                case 'M'
                    coneFundamentalsOptimizedForStimPosition = queryUserWhetherToEmployOptimizedConeFundamentals();
                    stimulusChromaticity = 'Mcone isolating';
                    responsesFileName = updateResponsesFileName(stimulusChromaticity, coneFundamentalsOptimizedForStimPosition, responsesFileName, identifierString);
                case 'S'
                    coneFundamentalsOptimizedForStimPosition = queryUserWhetherToEmployOptimizedConeFundamentals();
                    stimulusChromaticity = 'Scone isolating';
                    responsesFileName = updateResponsesFileName(stimulusChromaticity, coneFundamentalsOptimizedForStimPosition, responsesFileName, identifierString);
           end % switch
    end % while
end

function coneFundamentalsOptimizedForStimPosition = queryUserWhetherToEmployOptimizedConeFundamentals()
    coneFundamentalsOptimizedForStimPosition = 'invalid';
    
    while (~ismember(coneFundamentalsOptimizedForStimPosition, {'y', 'n'})) 
           coneFundamentalsOptimizedForStimPosition = lower(input('Employ cone fundamentals optimized for stimulus position? [y/n]: ', 's'));
    end

    coneFundamentalsOptimizedForStimPosition = strcmp(coneFundamentalsOptimizedForStimPosition, 'y');
    if (coneFundamentalsOptimizedForStimPosition)
        fprintf(2,'Employing cone fundamentals derived by probing the cone mosaic with monochromatic, equal energy stimuli delivered at the stimulus position\n');
    else
        fprintf(2,'Employing the default Stockman-Sharpe 2 deg cone fundamentals\n');
    end

end


function responsesFileName = updateResponsesFileName(stimulusChromaticity, coneFundamentalsOptimizedForStimPosition, responsesFileName, identifierString)
    
    if (coneFundamentalsOptimizedForStimPosition)
        coneFundamentalsString = 'optimized';
    else
        coneFundamentalsString = '';
    end

    switch (stimulusChromaticity)
        case 'achromatic'
            % We do not add a chromaticity in the filename if achromatic stimulus was employed - this is the default
        case 'Lcone isolating'
            responsesFileName = strrep(responsesFileName, identifierString, sprintf('%sLconeIsolating%s', coneFundamentalsString, identifierString));
        case 'Mcone isolating'
            responsesFileName = strrep(responsesFileName, identifierString, sprintf('%sMconeIsolating%s', coneFundamentalsString, identifierString));
        case 'Scone isolating'
            responsesFileName = strrep(responsesFileName, identifierString, sprintf('%sSconeIsolating%s', coneFundamentalsString, identifierString));
        otherwise
            error('Unknown stimulus chromaticity: ''%s''.', stimulusChromaticity);

     end % switch

end
