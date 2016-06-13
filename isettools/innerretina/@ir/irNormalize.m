function innerRetinaUn = irNormalize(innerRetinaUn, innerRetinaRef)

% mosaicNormalize

% Normalizes the linear response of one irPhys mosaic to another.


rLinSU = mosaicGet(innerRetinaUn.mosaic{1},'responseLinear');
rLin = mosaicGet(innerRetinaRef.mosaic{1},'responseLinear');
for cellNum = 1:length(rLinSU)
    maxTemp = max(rLin{cellNum,1} - innerRetinaRef.mosaic{1}.tonicDrive{cellNum});
    rLinSUNew{cellNum,1,1} = maxTemp*(rLinSU{cellNum,1} - innerRetinaUn.mosaic{1}.tonicDrive{cellNum})./max(rLinSU{cellNum,1}) + innerRetinaUn.mosaic{1}.tonicDrive{cellNum};% + 2.1;
%     rLinSUNew{cellNum,1,1} = [rLinSUNew{cellNum,1,1} min(rLin{cellNum,1})];
end
innerRetinaUn.mosaic{1} = mosaicSet(innerRetinaUn.mosaic{1},'responseLinear', rLinSUNew);