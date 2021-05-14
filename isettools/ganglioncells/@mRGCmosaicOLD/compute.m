% History:
%   11/06/20  dhb  Added noiseFactor key/value pair.

function  [mRGCresponses, temporalSupport] = compute(obj, coneMosaicResponses, timeAxis, varargin)
        p = inputParser;
        p.addParameter('seed', [], @isnumeric);
        p.addParameter('noiseFactor', [],@isnumeric);
        p.parse(varargin{:});
        
        % Apply the seed
        seed = p.Results.seed;
        if (~isempty(seed)), rng(seed); end
        
        % Set the noise factor
        if (~isempty(p.Results.noiseFactor))
            noiseFactor = p.Results.noiseFactor;
        else
            noiseFactor = obj.noiseFactor;
        end
        
        % Retrieve dimensions
        instancesNum = size(coneMosaicResponses,1);
        conesNum = size(coneMosaicResponses,2);
        rgcsNum = size(obj.coneWeights.center,2);
        if (ndims(coneMosaicResponses) == 2)
            nTimeBins = 1;
        elseif (ndims(coneMosaicResponses) == 3)
            nTimeBins = size(coneMosaicResponses,3);
        end
        
        % Ensure that the mosaic response dimensionality matches the
        % dimensionality of the input cone mosaic
        assert(conesNum == size(obj.inputConeMosaicMetaData.conePositionsDegs,1), ...
            sprintf('Dimensionality mismatch between cone mosaic and cone mosaic responses'));
        
        % Place holders for center/surround temporal impulse response functions
        centerImpulseResponse = [1];
        surroundImpulseResponse = [1];
        
        % Allocate memory
        mRGCresponses = zeros(instancesNum, rgcsNum, nTimeBins);
        temporalSupport = timeAxis;
        
        % Compute mRGC mosaic response
        for instanceIndex = 1:instancesNum
            
            % Extract response instance
            inputResponseInstance = squeeze(coneMosaicResponses(instanceIndex,:,:));
            if (nTimeBins == 1)
                inputResponseInstance = inputResponseInstance(:);
            end
            
            for RGCindex = 1:rgcsNum
                % Spatial integration
                centerConeWeights = (full(squeeze(obj.coneWeights.center(:, RGCindex))))';
                surroundConeWeights = (full(squeeze(obj.coneWeights.surround(:, RGCindex))))';
                centerPoolingResponse = centerConeWeights * inputResponseInstance;
                surroundPoolingResponse = surroundConeWeights * inputResponseInstance;
                
                % Temporal filtering
                centerResponse = conv(centerPoolingResponse, centerImpulseResponse, 'same');
                surroundResponse = conv(surroundPoolingResponse, surroundImpulseResponse, 'same');
                
                % Center-Surround opponency
                mRGCresponses(instanceIndex,RGCindex,:) = centerResponse - surroundResponse;
            end % RGCindex
            
        end % instanceIndex 
        
       
        % Noise
        if (~(strcmp(obj.noiseFlag, 'none'))) && (instancesNum>1)
            %fprintf(2,'Adding post-summation mRGC response noise (noiseFlag =''%s'').\n', obj.noiseFlag);
            % Mean over instances
            meanResponses = mean(mRGCresponses,1);
            % Max responses over time
            maxResponses = squeeze(max(meanResponses,3));
            
            % Add Gaussian noise with sigma = noiseFactor * max response
            for RGCindex = 1:rgcsNum
                sigma = maxResponses(RGCindex)*noiseFactor;
                mRGCresponses(:, RGCindex,:) = ...
                    meanResponses(:, RGCindex,:) + ...
                    randn(instancesNum,1,nTimeBins) * sigma;
            end
        else
            %fprintf(2,'No post-summation mRGC response noise (noiseFlag =''%s'').', obj.noiseFlag);
        end
        
end

