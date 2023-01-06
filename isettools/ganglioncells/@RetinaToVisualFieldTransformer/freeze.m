function freeze(obj)

    obj.theConeMosaic = [];
    obj.rfComputeStruct.modelConstants = [];

    obj.rfComputeStruct.theRotatedFittedVisualRF = [];
    obj.rfComputeStruct.theFittedVisualRFMap = [];
    obj.rfComputeStruct.theFittedVisualRFcenterConeMap = [];
    obj.rfComputeStruct.theFittedVisualRFsurroundConeMap = [];
    obj.rfComputeStruct.theRetinalRFcenterConeMap = [];
    obj.rfComputeStruct.theRetinalRFsurroundConeMap = [];
    obj.rfComputeStruct.targetVisualRFMap = [];
    obj.rfComputeStruct.targetVisualRFcenterMap = [];
    obj.rfComputeStruct.targetVisualRFsurroundMap = [];

end
