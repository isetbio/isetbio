function [noisyL noisyM noisyS] = coneSignalNoisyRepeats(ConeCurrentSignal, reps_max, coneSamplingRate, cone_mosaic)


[sz1 sz2 sz3] = size(ConeCurrentSignal);

ConeCurrentSignal_rs = reshape(ConeCurrentSignal, sz1*sz2, sz3);

for noise_reps = 1:reps_max
    ConeSignalFinalPlusNoise_rs = riekeAddNoise(ConeCurrentSignal_rs);%*coneSamplingRate);
    
    for cone_type=2:4
        
%         cone_locs = find(cone_mosaic==cone_type);
%         ConeSignalFinalCellPlusNoise{cone_type-1} = ConeSignalFinalPlusNoise_rs(cone_locs,:);
        
        % Only place the output signals corresponding to pixels in the mosaic
        % into the final output matrix.
        cone_locations = find(cone_mosaic==cone_type);
%         ConeSignal_rs = reshape(ConeSignalFinalPlusNoise_rs(:,:,1:sz3),[sz1*sz2],nSteps);
        ConeSignalFinal_rs(cone_locations,:) = ConeSignalFinalPlusNoise_rs(cone_locations,:);
        % obj.ConeSignalFinalCell{cone_type-1} = ConeSignal_rs(cone_locs,:);
        ConeSignalFinalCellPlusNoise{cone_type-1} = ConeSignalFinalPlusNoise_rs(cone_locations,:);
        
    end
    
%     % Reshape the output signal matrix.
%     obj.ConeCurrentSignal = reshape(ConeSignalFinal_rs, size(sensor.data.volts));
    
    
%     noisyL(noise_reps,:) = max(ConeSignalFinalCellPlusNoise{1});
%     noisyM(noise_reps,:) = max(ConeSignalFinalCellPlusNoise{2});
%     noisyS(noise_reps,:) = max(ConeSignalFinalCellPlusNoise{3});
    
    
    noisyL(noise_reps,:) = (ConeSignalFinalCellPlusNoise{1}(:));
    noisyM(noise_reps,:) = (ConeSignalFinalCellPlusNoise{2}(:));
    noisyS(noise_reps,:) = (ConeSignalFinalCellPlusNoise{3}(:));
    
    
    
end

