function activation2D = reshapeActivationMap1DTo2D(activation1D, pattern)

    nonNullCones = pattern(find(pattern>1));
    
    if (numel(activation1D) == numel(nonNullCones))
        iLsource = find(nonNullCones==2);
        iMsource = find(nonNullCones==3);
        iSsource = find(nonNullCones==4);

        iLdest = find(pattern==2);
        iMdest = find(pattern==3);
        iSdest = find(pattern==4);

        activation2D = zeros(size(pattern));
        activation2D(iLdest) = activation1D(iLsource);
        activation2D(iMdest) = activation1D(iMsource);
        activation2D(iSdest) = activation1D(iSsource);
     else
        error('The elements of the activation pattern (%d) does not match the number of non-null cones (%d) ', numel(activation), numel(nonNullCones));
    end
end
