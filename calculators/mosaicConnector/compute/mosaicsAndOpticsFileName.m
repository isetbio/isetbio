function [mosaicsFilename, opticsFilename, opticsPostFix,  PolansSubjectID] = mosaicsAndOpticsFileName(runParams)

    mosaicsFilename = sprintf('MosaicsForEccentricity_%2.0f_%2.0f_%2.0f_%2.0f_microns_coneSpecificity_%2.0f_orphanPolicy_%s.mat', ...
        runParams.rgcMosaicPatchEccMicrons(1), runParams.rgcMosaicPatchEccMicrons(2), ...
        runParams.rgcMosaicPatchSizeMicrons(1), runParams.rgcMosaicPatchSizeMicrons(2), ...
        runParams.maximizeConeSpecificity, runParams.orphanRGCpolicy);
    
    if ((runParams.noLCA == false) && (runParams.noOptics == false))
        opticsPostFix = 'normalOptics';
    elseif ((runParams.noLCA == true) && (runParams.noOptics == false))
        opticsPostFix = 'noLCAOptics';
    else
        opticsPostFix = 'noBlurOptics';
    end
    
    
    % Extract the Polans subjectID for the optics
    PolansSubjectID = runParams.deconvolutionOpticsParams.PolansWavefrontAberrationSubjectIDsToCompute;
    if (numel(PolansSubjectID )>1)
        PolansSubjectID = PolansSubjectID(1);
        fprintf(2, 'deconvolutionOpticsParams indicates more than 1 subject being used to derive the deconvolution model. Will generate optics for the first of these subjects\n');
    end
    
    opticsFilename = sprintf('OpticsForEccentricity_%2.0f_%2.0f_%2.0f_%2.0f_microns_PolansSID_%d_%s.mat', ...
        runParams.rgcMosaicPatchEccMicrons(1), runParams.rgcMosaicPatchEccMicrons(2), ...
        runParams.rgcMosaicPatchSizeMicrons(1), runParams.rgcMosaicPatchSizeMicrons(2), ...
        PolansSubjectID, opticsPostFix);
end
