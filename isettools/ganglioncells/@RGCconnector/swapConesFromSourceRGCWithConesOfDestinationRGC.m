function swapConesFromSourceRGCWithConesOfDestinationRGC(obj, ...
    sourceRGCconeIndicesToBeSwapped, sourceRGCindex, ...
    destinationRGCconeIndicesToBeSwapped, destinationRGCindex)


    for iSourceConeIndex = 1:numel(sourceRGCconeIndicesToBeSwapped)
        obj.transferConeFromSourceRGCToDestinationRGC( ...
            sourceRGCconeIndicesToBeSwapped(iSourceConeIndex), ...
            sourceRGCindex, destinationRGCindex);
    end

    for iDestinationConeIndex = 1:numel(destinationRGCconeIndicesToBeSwapped)
        obj.transferConeFromSourceRGCToDestinationRGC( ...
            destinationRGCconeIndicesToBeSwapped(iDestinationConeIndex), ...
            destinationRGCindex, sourceRGCindex);
    end

    
end

