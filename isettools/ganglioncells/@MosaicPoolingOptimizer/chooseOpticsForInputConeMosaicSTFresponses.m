function [opticsParams, opticsToEmploy, coneMosaicSTFresponsesFileName] = ...
    chooseOpticsForInputConeMosaicSTFresponses(mosaicParams)

    % Select optics to employ
    validOpticsChoices = {'n', 'c'};
    opticsChoice = 'invalid';
    while (~ismember(opticsChoice, validOpticsChoices))
       fprintf('\nOptics options: \n');
       fprintf('[n] - Native optics (at the mosaic''s center) \n');
       fprintf('[c] - Custom optics (at an arbitrary position) \n');
       opticsChoice = input('Optics to employ: ', 's');
    end

    % Generate the mosaic filename
    [mosaicFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('mosaic', ...
            'mosaicParams', mosaicParams);

    % Load the generated center-only connected mRGCmosaic
    load(fullfile(resourcesDirectory,mosaicFileName), 'theMidgetRGCMosaic');
    
    % Get the default optics params
    opticsParams = theMidgetRGCMosaic.defaultOpticsParams;

    % Generate filename for the computed coneMosaicSTF responses
    [coneMosaicSTFresponsesFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('coneMosaicSTFresponses', ...
            'mosaicParams', mosaicParams, ...
            'opticsParams', opticsParams);
   
    switch (opticsChoice)
        case 'n'
            opticsToEmploy = 'native';

        case 'c'
            opticsToEmploy = 'custom';

            % Get the custom optics position
            opticsParams.positionDegs = input('Enter (x,y) position (in degrees) for the optics (e.g., [2 0]) : ');
            opticsPositionPostfix = sprintf('AtCustomXYposition_%2.2f_%2.2f.mat', ...
                opticsParams.positionDegs(1), opticsParams.positionDegs(2));

            % Update the coneMosaicSTFresponsesFileName to reflect the custom optics position
            coneMosaicSTFresponsesFileName = strrep(coneMosaicSTFresponsesFileName, '.mat', opticsPositionPostfix);
            
        otherwise
            fprintf('Valid optics options: ''n'' (native, at mosaic''s center), or ''c'', (custom position)\n');
            error('Unknown optics choice: ''%s''.', opticsChoice)
    end

    coneMosaicSTFresponsesFileName = fullfile(resourcesDirectory, coneMosaicSTFresponsesFileName);
end
