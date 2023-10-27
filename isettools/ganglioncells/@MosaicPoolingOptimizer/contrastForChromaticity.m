function [coneContrasts, contrast] = contrastForChromaticity(stimulusChromaticity)

    switch (stimulusChromaticity)
        case 'achromatic'
            coneContrasts = [1 1 1];
            contrast = 0.75;
        case 'Lcone isolating'
            coneContrasts = [1 0 0];
            contrast = 0.12;
        case 'Mcone isolating'
            coneContrasts = [0 1 0];
            contrast = 0.12;
       case 'Scone isolating'
            coneContrasts = [0 0 1];
            contrast = 0.5;
        otherwise
            error('Unknown stimulus chromaticity: ''%s''.', stimulusChromaticity);
    end


end
