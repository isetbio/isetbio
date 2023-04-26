function arbitraryNodesToCompute = gridNodesToOptimize()
    arbitraryNodesToCompute = {};

    s = input('Compute all nodes (1) or speficic ones (2):');
    if (s == 1)
        return;
    else

        % Bad nodes
        for i = 15:29
        
         arbitraryNodesToCompute{numel(arbitraryNodesToCompute)+1} = struct(...
                'number', i, ...
                'coneType', [cMosaic.LCONE_ID cMosaic.MCONE_ID]);
        end


    end
end