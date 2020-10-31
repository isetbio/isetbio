function thePhotoCurrents = computeEccBasedPhotocurrentsFromIsomerizations(theConeMosaic, theConeMosaicMetaData, theIsomerizations, theIsomerizationsNull, rgcMosaicPatchEccMicrons)

    % Retrieve original data dimensions
    nTrials = size(theIsomerizations,1);
    conesNum = size(theIsomerizations,2);
    nDataPoints = size(theIsomerizations,3);
    
    % Save the full pattern
    theFullPattern = theConeMosaic.pattern;
    
    % Save the noise flag
    theNoiseFlag = theConeMosaic.os.noiseFlag;
    
    % change the mosaic pattern
    theConeMosaic.pattern = theConeMosaicMetaData.coneTypes;
    
    % Set the os eccentricityDegs property to match the ecc of the mosaic
    radialEccMicrons = sqrt(sum(rgcMosaicPatchEccMicrons.^2));
    theConeMosaic.os.eccentricityDegs = WatsonRGCModel.rhoMMsToDegs(radialEccMicrons/1000);
    theConeMosaic.os.noiseFlag = 'none';
    
    % Replicate the isomerizations sequence so that we allow enough time for
    % the photocurrent response to stabilize
    theIsomerizationsNull = cat(3, theIsomerizationsNull, theIsomerizationsNull);
    
    % Compute the mean current and the photocurrent IRs based on the null isomerizations
    [LMSfilters, meanCur] = theConeMosaic.computeCurrent(...
          'absorptionsInXWFormat', squeeze(theIsomerizationsNull(1,:,:)));
    
    theNullNoiseFreePhotocurrentResponse = theConeMosaic.current(:,end-nDataPoints+1:end);
    
    % Restore noise flag
    theConeMosaic.os.noiseFlag = theNoiseFlag;
    
    % Replicate the isomerizations sequence so that we allow enough time for
    % the photocurrent response to stabilize
    theIsomerizations = cat(3, theIsomerizations, theIsomerizations);
    
    % Compute photocurrent responses
    thePhotoCurrents = zeros(nTrials, conesNum, nDataPoints);
    parfor k = 1:nTrials
         theConeMosaic.computeCurrent(...
                    'interpFilters', LMSfilters, ...
                    'meanCur', meanCur, ...
                    'absorptionsInXWFormat', squeeze(theIsomerizations(k,:,:)));
         % Only keep the second half of the photocurrent responses
         % Return the modulation relative to the noise-free response to the NULL stimulus
         thePhotoCurrents(k,:,:) = bsxfun(@minus, theConeMosaic.current(:,end-nDataPoints+1:end), ...
                                   theNullNoiseFreePhotocurrentResponse(:,end));
    end
    
        
    % Restore the full pattern
    theConeMosaic.pattern = theFullPattern;
end

