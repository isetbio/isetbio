function retinalRFmodelParams = getSurroundParams(mosaicParams, opticsParams)

    retinalRFmodelParams = MosaicPoolingOptimizer.chooseRFmodelForSurroundConePoolingOptimization(...
        mosaicParams, opticsParams, ...
        'PackerDaceyH1horizontalCellIndex', 4);

end

