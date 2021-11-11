function electricallyCoupledResponses = electricallyCoupleConeResponses(obj, responses)

    % Compute neighoring cone indices and coupling weights with those
    % neighbors for each cone.
    
    if (obj.coneCouplingLambda < 0)
        % the 6 closest neigbors, but will only keep the closest one that is not an S-cone 
        maxCoupledNeigborsNum = 1+1;
    else
        % the 4 closest neighbors
        % Each cone contacts 4+/- 1.1 neighboring cones in the macaque fovea:
        % "Gap junctions between the pedicles of macaque foveal cones", 
        % Y. Tsukamoto, P. Masarachia, S. J. Schein, and P. Sterling
        % Vision Res. 32, 1809–1815 (1992).
        maxCoupledNeigborsNum = 1+4;
    end
    
    % For each cone, compute the coupling weights to its neighboring cones
    % and the indices of those cones only if these are empty 
    
    lookupNeighborsNum = 10;
    
    if (isempty(obj.coneCouplingWeights))
        
        % Compute and save indices of neighboring cones
        [~, neighboringConeIndices] = pdist2(obj.coneRFpositionsMicrons, obj.coneRFpositionsMicrons, ...
            'euclidean','Smallest', lookupNeighborsNum);
        
        obj_neighboringCoupledConeIndices = uint32(neighboringConeIndices');
        obj_coneRFpositionsMicrons = obj.coneRFpositionsMicrons;
        obj_coneCouplingLambda = obj.coneCouplingLambda;
        obj_coneApertureDiametersMicrons = obj.coneApertureDiametersMicrons;
        obj_coneTypes = obj.coneTypes;
        obj_coneCouplingWeights = zeros(size(obj.coneRFpositionsMicrons,1),lookupNeighborsNum, 'single');
        
        % Compute coupling weights
        for iCone = 1:size(obj.coneRFpositionsMicrons,1)
            
            neighboringConeIndices = obj_neighboringCoupledConeIndices(iCone,:);
            neighboringConePositionsMicrons = obj_coneRFpositionsMicrons(neighboringConeIndices,:);
                
            % Centered on the cone
            weightingFunctionCenterMicrons = obj_coneRFpositionsMicrons(iCone,:);
                
            % Weighting function
            neigboringConeDistancesMicrons = sqrt(sum((neighboringConePositionsMicrons - weightingFunctionCenterMicrons).^2,2));
            weights = exp(-(neigboringConeDistancesMicrons/(obj_coneCouplingLambda * obj_coneApertureDiametersMicrons(iCone))));
            
            
            % Weights for neighbors that are S-cones are 0. Coupling only between L/M-cones
            idx = find(obj_coneTypes(neighboringConeIndices) == cMosaic.SCONE_ID);
            weights(idx) = 0; 
                
 
            % If the target cone is an S-cone, all weights for neighboring cones are 0. Coupling only between L/M-cones
            if (obj_coneTypes(iCone) == cMosaic.SCONE_ID)
                 weights = weights * 0;
            end
            
            % Weight to the cone itself always 1.
            weights(1) = 1;
           
            % Save weights to the cone mosaic object
            obj_coneCouplingWeights(iCone,:) = single(weights);
        end
        
        obj.coneCouplingWeights  = obj_coneCouplingWeights;
        obj.neighboringCoupledConeIndices = obj_neighboringCoupledConeIndices;
        
        % Make coupling symmetrical
        for iCone = 1:size(obj.neighboringCoupledConeIndices,1)
            for iCoupledConeIndex = 2:size(obj.neighboringCoupledConeIndices,2)
                iCoupledCone = obj.neighboringCoupledConeIndices(iCone,iCoupledConeIndex);
                
                w1 = obj.coneCouplingWeights(iCone, iCoupledConeIndex);
                idx = find(obj.neighboringCoupledConeIndices(iCoupledCone,:) == iCone);
                w2 = obj.coneCouplingWeights(iCoupledCone, idx);
                
                if (obj.coneTypes(iCoupledCone) == cMosaic.SCONE_ID)
                    w1 = 0;
                    w2 = 0;
                end
                
                if (~isempty(w2))
                    wAverage = 0.5*(w1+w2);
                    obj_coneCouplingWeights(iCone, iCoupledConeIndex) = wAverage;
                    obj_coneCouplingWeights(iCoupledCone, idx) = wAverage;
                end
            end
        end
        

        obj.coneCouplingWeights = obj_coneCouplingWeights(:,1:maxCoupledNeigborsNum);
        obj.neighboringCoupledConeIndices = obj.neighboringCoupledConeIndices(:,1:maxCoupledNeigborsNum);
        
        if (obj.coneCouplingLambda < 0)
            obj.coneCouplingWeights = repmat([1 abs(obj.coneCouplingLambda)], [numel(obj_coneTypes) 1]);
            idx = find(obj_coneTypes == cMosaic.SCONE_ID);
            obj.coneCouplingWeights(idx,2) = 0;
        end
    end
    
    if (maxCoupledNeigborsNum == 2)
        % Single-shot coupling
        propagationStagesNum = 3;
    else
        % Propagate coupling in neighbors, mimicking syncytium
        propagationStagesNum = 3;
    end
    
    electricallyCoupledResponses = 0*responses;
    for stage = 1:propagationStagesNum
        
        fprintf('Coupling cone responses: stage %d of %d\n', stage, propagationStagesNum);
        for iCone = 1:size(obj.coneRFpositionsMicrons,1)
            neighboringConeIndices = obj.neighboringCoupledConeIndices(iCone,:);
            weights = obj.coneCouplingWeights(iCone,:);
            weights = reshape(weights, [1 1 numel(weights)]);
            electricallyCoupledResponses(:,:,iCone) = sum(bsxfun(@times, responses(:,:,neighboringConeIndices), weights),3);
        end
        responses = electricallyCoupledResponses;
    end
    
    
end




function [apertureKernel, theoreticalAreaMetersSquared, theoreticalAreaMicronsSquared, actualAreaMicronsSquared] = ...
                    coupledGaussianKernel(apertureRows, apertureCols, gaussianSigma, couplingLambda)
                  
    apertureRadii = sqrt(apertureRows .^ 2 + apertureCols .^ 2);
    
    % Gaussian kernel
    apertureKernel = exp(-0.5*(apertureRadii/gaussianSigma).^2);
    
    % Coupling kernel
    couplingKernel = exp(-(apertureRadii/couplingLambda));
    
    % Cascaded kernel
    apertureKernel = conv2(apertureKernel, couplingKernel, 'same');
    
    % Unit volume
    actualAreaMicronsSquared = sum(apertureKernel(:));
    apertureKernel = apertureKernel ./ actualAreaMicronsSquared;
    
    % Theoretical area based on continous space
    characteristicRadiusMicrons = gaussianSigma * sqrt(2.0);
    theoreticalAreaMetersSquared = ((pi * (characteristicRadiusMicrons*1e-6).^2));
    theoreticalAreaMicronsSquared = theoreticalAreaMetersSquared * 1e12;
    
    %Actual area based on discritized space
    dx = apertureCols(2)-apertureCols(1);
    actualAreaMicronsSquared = actualAreaMicronsSquared * dx^2;
end
