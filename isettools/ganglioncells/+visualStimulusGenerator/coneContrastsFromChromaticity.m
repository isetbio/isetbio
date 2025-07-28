function [coneContrasts, contrast] = coneContrastsFromChromaticity(stimulusChromaticity)

    switch (stimulusChromaticity)
        case 'Achromatic'
            coneContrasts = [1 1 1];
            contrast = 0.75;
        case 'LconeIsolating'
            coneContrasts = [1 0 0];
            contrast = 0.12;
        case 'MconeIsolating'
            coneContrasts = [0 1 0];
            contrast = 0.12;
       case 'SconeIsolating'
            coneContrasts = [0 0 1];
            contrast = 0.5;
        otherwise
            error('Unknown stimulus chromaticity: ''%s''.', stimulusChromaticity);
    end


end