function [sourceRGCinputs, sourceRGCweights, destRGCinputs, destRGCweights, ...
          sourceRGCcentroid, sourceRGCconeInputsNum, ...
          destRGCcentroid, destRGCconeInputsNum] = ...
                transferConeFromSourceRGCtoDestinationRGC(...
                        sourceRGCinputs, sourceRGCweights, ...
                        sourceRGCconeIndex, sourceRGCconeWeight, ...
                        destRGCinputs, destRGCweights, allConePositions)

    % Ensure that there are at least 2 cone inputs in the sourceRGC
    assert(numel(sourceRGCinputs) >= 2,'Source RGC must have at least 2 cone inputs.');
    
    % Remove sourceRGCconeIndex from the sourceRGC
    oldConeInputs = sourceRGCinputs;
    oldConeWeights = sourceRGCweights;
    [~, idx] = setdiff(oldConeInputs, sourceRGCconeIndex);
    oldConeInputs = oldConeInputs(idx);
    oldConeWeights = oldConeWeights(idx);
    sourceRGCinputs = oldConeInputs;
    sourceRGCweights = oldConeWeights;
  
    % Add sourceRGCconeIndex to the destinationRGC
    destRGCinputs(numel(destRGCinputs)+1) = sourceRGCconeIndex;
    destRGCweights(numel(destRGCweights)+1) = sourceRGCconeWeight;
    
    % Recompute centroids and cone inputs num to the source RGC
    [sourceRGCcentroid, sourceRGCconeInputsNum] = RGCRFconnector.centroidsFromConeInputs(...
                    {sourceRGCinputs}, {sourceRGCweights}, allConePositions);
                
    % Recompute centroids and cone inputs num to the destination RGC
    [destRGCcentroid, destRGCconeInputsNum] = RGCRFconnector.centroidsFromConeInputs(...
                    {destRGCinputs}, {destRGCweights}, allConePositions);

end