function [theOI, thePSF, theZcoeffs] = nativeOptics(obj, varargin)
    % Parse input
    p = inputParser;
    p.addParameter('atPosition', 'mosaic center', @(x)((ischar(x)&&(ismember(x, {'mosaic center', 'nearest RTVF object'})) || ((isnumeric(x))&&(numel(x)==2))) ));
    p.parse(varargin{:});

    if (ischar(p.Results.atPosition))
        switch (p.Results.atPosition) 
            case  'mosaic center'
                opticsPosition = obj.eccentricityDegs;
                d = bsxfun(@minus,obj.multifocalRTVFgrids.samplingPositionGrid, obj.eccentricityDegs);
                [~, iObj] = min(sum(d.^2,2));
                opticsParams = obj.multifocalRTVFopticsParams{iObj};

            case 'nearest RTVF object'
                d = bsxfun(@minus,obj.multifocalRTVFgrids.samplingPositionGrid, obj.eccentricityDegs);
                [~, iObj] = min(sum(d.^2,2));
                opticsParams = obj.multifocalRTVFopticsParams{iObj};
                opticsPosition = obj.multifocalRTVFgrids.samplingPositionGrid(iObj,:);

            otherwise
                error('Unknown optics position: ''%s''.', p.Results.atPosition);
        end
    else
        opticsPosition = p.Results.atPosition;
        d = bsxfun(@minus,obj.multifocalRTVFgrids.samplingPositionGrid, opticsPosition);
        [~, iObj] = min(sum(d.^2,2));
        opticsParams = obj.multifocalRTVFopticsParams{iObj};
    end

    % Generate the OI based on the retrieved opticsParams
    [oiEnsemble, psfEnsemble, theZcoeffs] = obj.inputConeMosaic.oiEnsembleGenerate(opticsPosition, ...
                        'zernikeDataBase', opticsParams.ZernikeDataBase, ...
                        'subjectID', opticsParams.testSubjectID, ...
                        'pupilDiameterMM', opticsParams.pupilDiameterMM, ...
                        'refractiveErrorDiopters', opticsParams.refractiveErrorDiopters, ...
                        'zeroCenterPSF', opticsParams.zeroCenterPSF, ...
                        'subtractCentralRefraction', opticsParams.subtractCentralRefraction, ...
                        'wavefrontSpatialSamples', opticsParams.wavefrontSpatialSamples, ...
                        'upsampleFactor', opticsParams.psfUpsampleFactor, ...
                        'warningInsteadOfErrorForBadZernikeCoeffs', true);    

    % Return theOI and thePSF
    theOI = oiEnsemble{1};
    thePSF = psfEnsemble{1};
end