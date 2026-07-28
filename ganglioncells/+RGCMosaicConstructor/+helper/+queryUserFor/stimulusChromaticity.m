function [stimulusChromaticity, coneFundamentalsOptimizedForStimPosition, theUpdatedFileName] = ...
        stimulusChromaticity(theFilename, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('doNotQueryUserInsteadUseThisChromaParams', [], @(x)(isempty(x) || isstruct(x)));
    p.parse(varargin{:});
    chromaParamsToEmploy = p.Results.doNotQueryUserInsteadUseThisChromaParams;

    if (~isempty(chromaParamsToEmploy))
        % Used passed params
        stimulusChromaticity = chromaParamsToEmploy.stimulusChromaticity;
        coneFundamentalsOptimizedForStimPosition = chromaParamsToEmploy.coneFundamentalsOptimizedForStimPosition;   
    else
        % Query user
        chromaChoice = 'invalid';
        validChromaticityChoices = { ...
    		'A', ...
            'L', ...
            'M', ...
            'S' };

    	while (~ismember(upper(chromaChoice), validChromaticityChoices)) 
         	fprintf('\nChoose stimulus chromaticity \n');
            fprintf('[A]  - Achromatic (default) \n');
            fprintf('[L]  - L-cone isolating \n');
            fprintf('[M]  - M-cone isolating \n');
            fprintf('[S]  - S-cone isolating \n');
            chromaChoice = input('Stimulus chromaticity: ', 's');
       	end

       	switch (upper(chromaChoice))
            case 'A'
                stimulusChromaticity = 'Achromatic';
                coneFundamentalsOptimizedForStimPosition = queryUserWhetherToEmployOptimizedConeFundamentals();
            case 'L'
                coneFundamentalsOptimizedForStimPosition = queryUserWhetherToEmployOptimizedConeFundamentals();
                stimulusChromaticity = 'LconeIsolating'; 
            case 'M'
                coneFundamentalsOptimizedForStimPosition = queryUserWhetherToEmployOptimizedConeFundamentals();
                stimulusChromaticity = 'MconeIsolating';
            case 'S'
                coneFundamentalsOptimizedForStimPosition = queryUserWhetherToEmployOptimizedConeFundamentals();
                stimulusChromaticity = 'SconeIsolating';
        end % switch
    end

    if (coneFundamentalsOptimizedForStimPosition)
        coneFundamentalsString = sprintf('optimized_%s', stimulusChromaticity);
    else
        coneFundamentalsString = sprintf('defaultSS_%s', stimulusChromaticity);
    end

    theUpdatedFileName = strrep(theFilename, '.mat', sprintf('_%s.mat', coneFundamentalsString));
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