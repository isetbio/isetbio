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

    resourcesDirectory = fullfile(...
                MosaicPoolingOptimizer.localDropboxPath, ...
                'productionMidgetRGCMosaics/MosaicOptimizerResources');

    pdfsDirectory = fullfile(...
                MosaicPoolingOptimizer.localDropboxPath, ...
                'productionMidgetRGCMosaics/MosaicOptimizerPDFs');

    switch (component)
        case 'mosaic'
            resourceFileName = generateMosaicFileName(mosaicParams);
            
        case 'coneMosaicSTFresponses'
            resourceFileName = generateConeMosaicSTFResponsesFileName(mosaicParams, opticsParams);

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
    if (isempty(mosaicParams))
        error('Must provide mosaicParams struct')
    end
    m = sprintf('mRGCMosaicEcDegs(%2.1f_%2.1f)_SizeDegs(%2.1f_%2.1f).mat', ...
             mosaicParams.eccDegs(1), mosaicParams.eccDegs(2), ...
             mosaicParams.sizeDegs(1), mosaicParams.sizeDegs(2));
end

function m = generateConeMosaicSTFResponsesFileName(mosaicParams, opticsParams)
    mosaicFileName = generateMosaicFileName(mosaicParams);
    opticsString = generateOpticsString(opticsParams);
    m = strrep(mosaicFileName, '.mat', sprintf('%s_coneMosaicSTFresponses.mat', opticsString));
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
    m = sprintf('_%s_SubjRank%d_%s_%2.2fmmPupil_RefrErrorDiopters_%2.2f', ...
        opticsParams.ZernikeDataBase, ...
        opticsParams.examinedSubjectRankOrder, ...
        strrep(opticsParams.analyzedEye, ' ', '_'),...
        opticsParams.pupilDiameterMM, ...
        opticsParams.refractiveErrorDiopters);
end