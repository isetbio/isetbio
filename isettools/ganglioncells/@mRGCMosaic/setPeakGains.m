% Method to set the peak gain of each mRGC 
function setPeakGains(obj, method, methodParams)
    switch (method)
        case 'arbitrary'
            obj.rgcRFgains = methodParams;

        case 'CK95formulaAppliedToMRGCMosaicVisualRcDegs'
            % From Croner&Kaplan '95, Figure 5b
            % Here, the passed methodParams variable must be computed externally
            % as a [1 x rgcsNum] vector of the visual Rc (degs) obtained by 
            % fitting the visual STF
            visualRcDegs = methodParams;
            obj.rgcRFgains = 0.391 * visualRcDegs .^ (-1.85);

        case '1/integrated center cone weights'
            if (isempty(obj.rgcRFcenterConePoolingMatrix))
                error('Cannot call the setPeakGains method before the rgcRFcenterConePoolingMatrix is computed.');
            else
                centerConeWeightsSum = zeros(1,obj.rgcsNum);
                parfor iRGC = 1:obj.rgcsNum
                    connectivityVector = full(squeeze(obj.rgcRFcenterConePoolingMatrix(:, iRGC)));
                    centerConeIndices = find(abs(connectivityVector) > 0.0001);
                    centerConeWeightsSum(iRGC) = sum(connectivityVector(centerConeIndices));
                end

                maxGain = methodParams(1);
                obj.rgcRFgains = maxGain  * (centerConeWeightsSum .^ (-1.0));
            end

        case '1/integrated center retinal cone apertures'
            if (isempty(obj.rgcRFcenterConePoolingMatrix))
                error('Cannot call the setPeakGains method before the rgcRFcenterConePoolingMatrix is computed.');
            else
                integratedConeApertureAreas = zeros(1,obj.rgcsNum);
                parfor iRGC = 1:obj.rgcsNum
                    connectivityVector = full(squeeze(obj.rgcRFcenterConePoolingMatrix(:, iRGC)));
                    centerConeIndices = find(abs(connectivityVector) > 0.0001);
                    integratedConeApertureAreas(iRGC) = sum(obj.inputConeMosaic.computeApertureAreasMetersSquared(centerConeIndices));
                end
                
                maxGain = methodParams(1);
                obj.rgcRFgains = maxGain * (integratedConeApertureAreas .^ (-1.0));
            end

        otherwise
            fprintf(2, '\nValid peak gain setting methods: {''CK95formulaAppliedToMRGCMosaicVisualRcDegs'', ''1/integrated center cone weights'', ''1/integrated center retinal cone apertures''}\n')
            error('Invalid peak gain setting method: ''%s''.', method);
    end % switch

end