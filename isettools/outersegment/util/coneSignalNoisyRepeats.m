function [noisyL noisyM noisyS] = coneSignalNoisyRepeats(ConeCurrentSignal, reps_max, coneSamplingRate, cone_mosaic)


[sz1 sz2 sz3] = size(ConeCurrentSignal);

ConeCurrentSignal_rs = reshape(ConeCurrentSignal, sz1*sz2, sz3);

for noise_reps = 1:reps_max
    ConeSignalFinalPlusNoise_rs = riekeAddNoise(ConeCurrentSignal_rs);%*coneSamplingRate);
    
    for cone_type=2:4
        
        cone_locs = find(cone_mosaic==cone_type);
        ConeSignalFinalCellPlusNoise{cone_type-1} = ConeSignalFinalPlusNoise_rs(cone_locs,:);
    end
    
    noisyL(noise_reps,:) = max(ConeSignalFinalCellPlusNoise{1});
    noisyM(noise_reps,:) = max(ConeSignalFinalCellPlusNoise{2});
    noisyS(noise_reps,:) = max(ConeSignalFinalCellPlusNoise{3});
    
      
    
end

