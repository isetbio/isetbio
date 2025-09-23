function sampledCenterConesNum = downsampledRFcenterConeNumerosityPopulation(uniqueCenterConesNum)
    
    minConesNum = min(uniqueCenterConesNum);
    maxConesNum = max(uniqueCenterConesNum);


    downsampledConesNum = 1:9;
    downsampledConesNum = cat(2,downsampledConesNum, 11:2:19);
    downsampledConesNum = cat(2,downsampledConesNum, 21:3:30);
    downsampledConesNum = cat(2,downsampledConesNum, 34:4:50);
    downsampledConesNum = cat(2,downsampledConesNum, 55:5:70);
    downsampledConesNum = cat(2,downsampledConesNum, 80:10:100);
    downsampledConesNum = cat(2,downsampledConesNum, 120:20:200);
    downsampledConesNum = cat(2,downsampledConesNum, 230:30:320);
    downsampledConesNum = cat(2,downsampledConesNum, 360:40:600);
    downsampledConesNum = cat(2,downsampledConesNum, 650:50:2000);
    downsampledConesNum = cat(2,downsampledConesNum, 2100:100:4000);
    downsampledConesNum = cat(2,downsampledConesNum, 4200:200:10000);
    downsampledConesNum = cat(2,downsampledConesNum, 10500:500:maxConesNum);

    downsampledConesNum = downsampledConesNum(...
        (downsampledConesNum >= minConesNum) & ...
        (downsampledConesNum <= maxConesNum)...
        );
   
    uniqueCenterConesNumWithinMappedDownSampledCenterCones = cell(1, numel(downsampledConesNum));
    theMappedDownSampledCenterConesNum = [];

    % Map each of the uniqueCenterConesNum to one of the downsampledConesNum
    for i = 1:numel(uniqueCenterConesNum)
        theUniqueCenterConesNum = uniqueCenterConesNum(i);

        % Find theMappedSampledCenterConesNum for theUniqueCenterConesNum
        [d,idx] = min(abs(theUniqueCenterConesNum-downsampledConesNum));
        theMappedDownSampledCenterConesNum(numel(theMappedDownSampledCenterConesNum)+1) = downsampledConesNum(idx);

        % Collect all theUniqueCenterConesNum that are mapped to theMappedDownSampledCenterConesNum
        dd = uniqueCenterConesNumWithinMappedDownSampledCenterCones{idx};
        if (isempty(dd))
            dd = theUniqueCenterConesNum;
        else
            dd = cat(1,dd,theUniqueCenterConesNum);
        end
        uniqueCenterConesNumWithinMappedDownSampledCenterCones{idx} = dd;
    end

    
    sampledCenterConesNum = [];
    for i = 1:numel(uniqueCenterConesNumWithinMappedDownSampledCenterCones)
        theUniqueCenterConesNumsWithinThisDownsampledCenterConesNum = uniqueCenterConesNumWithinMappedDownSampledCenterCones{i};
        if (numel(theUniqueCenterConesNumsWithinThisDownsampledCenterConesNum) > 0)
            % Select the mean uniqueCenterConesNum 
            [~,idx] = min(abs(theUniqueCenterConesNumsWithinThisDownsampledCenterConesNum - mean(theUniqueCenterConesNumsWithinThisDownsampledCenterConesNum)));
            sampledCenterConesNum(numel(sampledCenterConesNum)+1) = theUniqueCenterConesNumsWithinThisDownsampledCenterConesNum(idx);
        end
    end

    if (1==2)
        demoMode = false;
        if (demoMode)
            uniqueCenterConesNum
            theMappedDownSampledCenterConesNum
            sampledCenterConesNum = unique(sampledCenterConesNum)
            pause
        else
            fprintf('\nUnique center cones num in this patch:\n');
            uniqueCenterConesNum(:)
            fprintf('\nOptimized center cones num:\n');
            sampledCenterConesNum(:)
        end
    end

end
