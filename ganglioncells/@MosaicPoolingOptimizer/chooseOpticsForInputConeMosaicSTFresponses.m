function [opticsParams, opticsToEmploy, coneMosaicSTFresponsesFileName] = ...
    chooseOpticsForInputConeMosaicSTFresponses(mosaicParams, varargin )

    validOpticsChoices = {'n', 'c', 'ao'};

    p = inputParser;
    p.addParameter('opticsChoice', [], @(x)(isempty(x)||ismember(x,validOpticsChoices)));
    p.addParameter('refractiveErrorDiopters', [], @isscalar);
    p.addParameter('queryUserWhetherToUseSpecializedOpticsCases', true, @islogical);

    p.parse(varargin{:});
    opticsChoice = p.Results.opticsChoice;
    refractiveErrorDiopters = p.Results.refractiveErrorDiopters;
    queryUserWhetherToUseSpecializedOpticsCases = p.Results.queryUserWhetherToUseSpecializedOpticsCases;

    % Select optics to employ
    if (isempty(opticsChoice))
        opticsChoice = 'invalid';
        while (~ismember(opticsChoice, validOpticsChoices)) 
           fprintf('\nOptics options for mosaic at (%2.1f,%2.1f): \n', mosaicParams.eccDegs(1), mosaicParams.eccDegs(2));
           fprintf('[n]  - Native optics (at the mosaic''s center) \n');
           fprintf('[c]  - Custom optics (at an arbitrary position) \n');
           fprintf('[ao] - Adaptive optics \n');
           opticsChoice = input('Optics to employ: ', 's');
        end
    end

    % Generate the mosaic filename
    [mosaicFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('mosaic', ...
            'mosaicParams', mosaicParams);

    % Load the generated center-only connected mRGCmosaic
    load(fullfile(resourcesDirectory,mosaicFileName), 'theMidgetRGCMosaic');
    
    % Get the default optics params
    opticsParams = theMidgetRGCMosaic.defaultOpticsParams;

    % Ask user whether to override the optics with a vLambda-weighted
    % monochromatic PSF
    opticsParams.employMonochromaticVlambdaWeightedPSF = false;

    if (queryUserWhetherToUseSpecializedOpticsCases)
        tmp = lower(input('Override optics with monochromatic, vLambda-weighted PSF? [y/n]: ', 's'));
        if (strcmp(tmp, 'y'))
            opticsParams.employMonochromaticVlambdaWeightedPSF = true;
        end
    end

    % Generate filename for the computed coneMosaicSTF responses
    [coneMosaicSTFresponsesFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('coneMosaicSTFresponses', ...
            'mosaicParams', mosaicParams, ...
            'opticsParams', opticsParams);
   
    switch (opticsChoice)
        case 'n'
            opticsToEmploy = 'native';

            if (isempty(refractiveErrorDiopters))
                opticsParams.refractiveErrorDiopters = input('Enter refractive error in diopters, e.g. 0.15 : ');
                if (isempty(opticsParams.refractiveErrorDiopters))
                    opticsParams.refractiveErrorDiopters = 0.0;
                end
            else
                opticsParams.refractiveErrorDiopters = refractiveErrorDiopters;
            end


            if (abs(opticsParams.refractiveErrorDiopters) > 0)
                % Get updated coneMosaicSTFresponsesFileName
                coneMosaicSTFresponsesFileName = MosaicPoolingOptimizer.resourceFileNameAndPath('coneMosaicSTFresponses', ...
                    'mosaicParams', mosaicParams, ...
                    'opticsParams', opticsParams);
            end

        case 'c'
            opticsToEmploy = 'custom';

            % Get the custom optics position
            opticsParams.positionDegs = input('Enter (x,y) position (in degrees) for the optics (e.g., [2 0]) : ');
            opticsPositionPostfix = sprintf('AtCustomXYposition_%2.2f_%2.2f.mat', ...
                opticsParams.positionDegs(1), opticsParams.positionDegs(2));

            % Update the coneMosaicSTFresponsesFileName to reflect the custom optics position
            coneMosaicSTFresponsesFileName = strrep(coneMosaicSTFresponsesFileName, '.mat', opticsPositionPostfix);
            
        case 'ao'
            opticsToEmploy = 'adaptive optics';

            % Set the subject to 0 to indicate AO
            opticsParams.examinedSubjectRankOrder = 0;

            % No LCA for AO
            opticsParams.noLCA = true;

            % High upsample factor to capture the diffraction-limited PSF
            opticsParams.psfUpsampleFactor = 3;
            opticsParams.wavefrontSpatialSamples = 301;

            % Query user for AO pupil diameter
            opticsParams.pupilDiameterMM = input('Enter AO pupil diameter in mm, e.g., 6 : ');
            if (isempty(opticsParams.pupilDiameterMM))
                opticsParams.pupilDiameterMM = 6.0;
            end
            assert((opticsParams.pupilDiameterMM > 0) && (opticsParams.pupilDiameterMM <=6), ...
                'pupil diameter must be greater than 0 and not larger than 6.0 mm');

            % Query user for refractive error 
            opticsParams.refractiveErrorDiopters = input('Enter refractive error in diopters, e.g. 0.15 : ');
            if (isempty(opticsParams.refractiveErrorDiopters))
                opticsParams.refractiveErrorDiopters = 0.0;
            end

            % Get updated coneMosaicSTFresponsesFileName
            coneMosaicSTFresponsesFileName = ...
                MosaicPoolingOptimizer.resourceFileNameAndPath('coneMosaicSTFresponses', ...
                'mosaicParams', mosaicParams, ...
                'opticsParams', opticsParams);

        otherwise
            fprintf('Valid optics options: ''n'' (native, at mosaic''s center), ''c'', (custom position), or ''ao'' (adaptive optics)\n');
            error('Unknown optics choice: ''%s''.', opticsChoice)
    end

    coneMosaicSTFresponsesFileName = fullfile(resourcesDirectory, coneMosaicSTFresponsesFileName);
end
