function centerConesNumCases = centerConePoolingStats(obj)
    
    if (isempty(obj.rgcRFcenterConePoolingMatrix))
        centerConesNum = full(sum(obj.rgcRFcenterConeConnectivityMatrix,1));
    else
        centerConesNum = full(sum(obj.rgcRFcenterConePoolingMatrix,1));
    end

    centerConesNumCases = unique(centerConesNum);
    fprintf('\n* * * mRGCmosaic stats * * * \n');
    for iCase = 1:numel(centerConesNumCases)
        thisManyRGCs = numel(find(centerConesNum == centerConesNumCases(iCase)));
        fprintf('There are %05d RGCs with %d-cone(s) in their RF center\n', ...
            thisManyRGCs, centerConesNumCases(iCase));
    end
    fprintf('\n* * * * * * * * * * * * * * * \n');

end

