function [resourceFileName, resourcesDirectory, pdfsDirectory] = resourceFileNameAndPath(component, varargin)
    % Parse input
    p = inputParser;
    p.addParameter('mosaicParams', [], @(x)(isempty(x)||(isstruct(x))));
    p.addParameter('opticsParams', [], @(x)(isempty(x)||(isstruct(x))));
    p.addParameter('retinalRFmodelParams', [], @(x)(isempty(x)||(isstruct(x))));

    p.parse(varargin{:});
    mosaicParams = p.Results.mosaicParams;
    opticsParams = p.Results.opticsParams;
    retinalRFmodelParams = p.Results.retinalRFmodelParams;
    

    % Validate mosaicParams.rgcType
    assert(ismember(mosaicParams.rgcType, {'ONcenterMidgetRGC'}), ...
        sprintf('Invalid mosaicParams.rgcType: ''%s''.', mosaicParams.rgcType));

    resourcesRootDirectory = fullfile(...
                MosaicPoolingOptimizer.localDropboxPath, ...
                sprintf('%smosaics',mosaicParams.rgcType));

    pdfsDirectory = fullfile(...
                resourcesRootDirectory, ...
                'pdfs');

    if (strcmp(component, 'computeReadyMosaic'))
        resourcesDirectory = fullfile(...
                resourcesRootDirectory, ...
                'computeReadyMosaics');
    else
        resourcesDirectory = fullfile(...
                resourcesRootDirectory, ...
                'intermediateFiles');
    end


    
    switch (component)
        case 'mosaic'
            resourceFileName = generateMosaicFileName(mosaicParams);
            
        case 'coneMosaicSTFresponses'
            resourceFileName = generateConeMosaicSTFResponsesFileName(mosaicParams, opticsParams);

        case 'mRGCMosaicSTFresponses'
            resourceFileName = generateMRGCMosaicSTFResponsesFileName(mosaicParams, opticsParams);

        case 'coneMosaicSubspaceResponses'
            resourceFileName = generateConeMosaicSubspaceResponsesFileName(mosaicParams, opticsParams);

        case 'mRGCMosaicSubspaceResponses'
            resourceFileName = generateMRGCMosaicSubspaceResponsesFileName(mosaicParams, opticsParams);

        case 'coneMosaicMSequenceResponses'
            resourceFileName = generateConeMosaicMSequenceResponsesFileName(mosaicParams, opticsParams);

        case 'mRGCMosaicMSequenceResponses'
            resourceFileName = generateMRGCMosaicMSequenceResponsesFileName(mosaicParams, opticsParams);

        case 'optimizedRGCpoolingObjects'
            resourceFileName = generateOptimizedRGCpoolingObjectsFileName(mosaicParams, opticsParams, retinalRFmodelParams);

        case 'computeReadyMosaic'
            resourceFileName = generateComputeReadyMidgetRGCMosaicFileName(mosaicParams, opticsParams, retinalRFmodelParams);

        case 'pdfsDirectory'
            resourceFileName = '';

        otherwise
            error('Unknown component: ''%s''.', component)
    end

end

function m = generateMosaicFileName(mosaicParams)
    % Validate mosaicParams
    validateMosaicParams(mosaicParams);

    m = sprintf('mRGCMosaicEcDegs(%2.1f_%2.1f)_SizeDegs(%2.1f_%2.1f).mat', ...
             mosaicParams.eccDegs(1), mosaicParams.eccDegs(2), ...
             mosaicParams.sizeDegs(1), mosaicParams.sizeDegs(2));
end

function validateMosaicParams(mosaicParams)
    if (isempty(mosaicParams))
        error('Must provide a non-empty mosaicParams struct with ''eccDegs'' and ''sizeDegs'' fields.')
    end
    
    assert(isfield(mosaicParams, 'eccDegs'), ...
        sprintf('''mosaicParams'': missing field ''eccDegs'', containing the [x,y] center coords for the mosaic.'));

    assert(numel(mosaicParams.eccDegs) == 2, ...
        sprintf('''mosaicParams.eccDegs'': must contain a 2-element vector with the [x,y] center coords for the mosaic.'));

    assert(isfield(mosaicParams, 'sizeDegs'), ...
        sprintf('''mosaicParams'': missing field ''sizeDegs'', containing the [w,h] size of the mosaic.'));

    assert(numel(mosaicParams.sizeDegs) == 2, ...
        sprintf('''mosaicParams.sizeDegs'': must contain a 2-element vector with the [w,h] size for the mosaic.'));

end

function m = generateConeMosaicSTFResponsesFileName(mosaicParams, opticsParams)
    mosaicFileName = generateMosaicFileName(mosaicParams);
    opticsString = generateOpticsString(opticsParams);
    m = strrep(mosaicFileName, '.mat', sprintf('%s_coneMosaicSTFresponses.mat', opticsString));
end

function m = generateMRGCMosaicSTFResponsesFileName(mosaicParams, opticsParams)
    mosaicFileName = generateMosaicFileName(mosaicParams);
    opticsString = generateOpticsString(opticsParams);
    m = strrep(mosaicFileName, '.mat', sprintf('%s_mRGCMosaicSTFresponses.mat', opticsString));
end

function m = generateConeMosaicSubspaceResponsesFileName(mosaicParams, opticsParams)
    mosaicFileName = generateMosaicFileName(mosaicParams);
    opticsString = generateOpticsString(opticsParams);
    m = strrep(mosaicFileName, '.mat', sprintf('%s_coneMosaicSubspaceResponses.mat', opticsString));
end

function m = generateMRGCMosaicSubspaceResponsesFileName(mosaicParams, opticsParams)
    mosaicFileName = generateMosaicFileName(mosaicParams);
    opticsString = generateOpticsString(opticsParams);
    m = strrep(mosaicFileName, '.mat', sprintf('%s_mRGCMosaicSubspaceResponses.mat', opticsString));
end

function m = generateConeMosaicMSequenceResponsesFileName(mosaicParams, opticsParams)
    mosaicFileName = generateMosaicFileName(mosaicParams);
    opticsString = generateOpticsString(opticsParams);
    m = strrep(mosaicFileName, '.mat', sprintf('%s_coneMosaicMSequenceResponses.mat', opticsString));
end

function m = generateMRGCMosaicMSequenceResponsesFileName(mosaicParams, opticsParams)
    mosaicFileName = generateMosaicFileName(mosaicParams);
    opticsString = generateOpticsString(opticsParams);
    m = strrep(mosaicFileName, '.mat', sprintf('%s_mRGCMosaicMSequenceResponses.mat', opticsString));
end

function m = generateOptimizedRGCpoolingObjectsFileName(mosaicParams, opticsParams, retinalRFmodelParams)
    mosaicFileName = generateMosaicFileName(mosaicParams);
    opticsString = generateOpticsString(opticsParams);
    if (contains(retinalRFmodelParams.conePoolingModel, 'H1cellIndex1'))
        H1cellIndex = 1;
    elseif (contains(retinalRFmodelParams.conePoolingModel, 'H1cellIndex2'))
        H1cellIndex = 2;
    elseif (contains(retinalRFmodelParams.conePoolingModel, 'H1cellIndex3'))
        H1cellIndex = 3;
    elseif (contains(retinalRFmodelParams.conePoolingModel, 'H1cellIndex4'))
        H1cellIndex = 4;
    else
        error('Could not determine H1cellIndex\n');
    end

    m = strrep(mosaicFileName, '.mat', sprintf('%s_H1cellIndex%dBasedConePoolingObject.mat', opticsString, H1cellIndex));
end


function m = generateComputeReadyMidgetRGCMosaicFileName(mosaicParams, opticsParams, retinalRFmodelParams)
    mosaicFileName = generateMosaicFileName(mosaicParams);
    opticsString = generateOpticsString(opticsParams);
    if (contains(retinalRFmodelParams.conePoolingModel, 'H1cellIndex1'))
        H1cellIndex = 1;
    elseif (contains(retinalRFmodelParams.conePoolingModel, 'H1cellIndex2'))
        H1cellIndex = 2;
    elseif (contains(retinalRFmodelParams.conePoolingModel, 'H1cellIndex3'))
        H1cellIndex = 3;
    elseif (contains(retinalRFmodelParams.conePoolingModel, 'H1cellIndex4'))
        H1cellIndex = 4;
    else
        error('Could not determine H1cellIndex\n');
    end

    m = strrep(mosaicFileName, '.mat', sprintf('%s_H1cellIndex%d_ComputeReadyMosaic.mat', opticsString, H1cellIndex));
end

function m = generateOpticsString(opticsParams)
    if (isempty(opticsParams.positionDegs))
        if (opticsParams.examinedSubjectRankOrder == 0)
            % AO optics
            m = sprintf('_adaptiveOptics_%2.2fmmPupil_ResidualDefocusDiopters_%2.2f', ...
                opticsParams.pupilDiameterMM, ...
                opticsParams.refractiveErrorDiopters);
        else
            % Human subject optics
            m = sprintf('_%s_SubjRank%d_%s_%2.2fmmPupil_RefrErrorDiopters_%2.2f', ...
                opticsParams.ZernikeDataBase, ...
                opticsParams.examinedSubjectRankOrder, ...
                strrep(opticsParams.analyzedEye, ' ', '_'),...
                opticsParams.pupilDiameterMM, ...
                opticsParams.refractiveErrorDiopters);
        end
    else
        % Human subject optics
        m = sprintf('_%s_SubjRank%d_%s_atCustomXYposition_%2.2f_%2.2f_%2.2fmmPupil_RefrErrorDiopters_%2.2f', ...
            opticsParams.ZernikeDataBase, ...
            opticsParams.examinedSubjectRankOrder, ...
            strrep(opticsParams.analyzedEye, ' ', '_'),...
            opticsParams.pupilDiameterMM, ...
            opticsParams.refractiveErrorDiopters);
    end

    if (isfield(opticsParams, 'employMonochromaticVlambdaWeightedPSF'))
        if (opticsParams.employMonochromaticVlambdaWeightedPSF)
            m = sprintf('%s_MonoVlambdaPSF', m);
        end
    end
end
