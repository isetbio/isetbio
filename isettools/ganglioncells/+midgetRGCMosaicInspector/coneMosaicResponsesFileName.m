function  responsesFileName = coneMosaicResponsesFileName(mosaicFileName, opticsPositionDegs)

    RGCMosaicResponsesFileName = midgetRGCMosaicInspector.responsesFileName(...
        mosaicFileName, opticsPositionDegs);
    
    responsesFileName = strrep(RGCMosaicResponsesFileName, 'Responses', 'InputConeMosaicResponses');
end
