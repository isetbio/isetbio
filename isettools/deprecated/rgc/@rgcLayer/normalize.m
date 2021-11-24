function innerRetinaUn = normalize(innerRetinaUn, innerRetinaRef)
% Normalize an inner retina object
%
% Syntax:
%    innerRetinaUn = normalize(innerRetinaUn, innerRetinaRef)
%
% Description:
%    Normalize the linear response of one irPhys mosaic to another.
%
% Inputs:
%    innerRetinaUn  - Object. An irPhys mosaic object.
%    innerRetinaRef - Object. A second irPhys mosaic object (to normalize
%                     innerRetinaUn to).
%
% Outputs:
%    innerRetinaUn  - Object. The modified irPhys mosaic object.
%
% Optional key/value pairs:
%    None.
%

rLinSU = mosaicGet(innerRetinaUn.mosaic{1}, 'responseLinear');
rLin = mosaicGet(innerRetinaRef.mosaic{1}, 'responseLinear');
for cellNum = 1:size(rLinSU, 1)
    rLinVar(cellNum, 1, :) = rLin(cellNum, 1, :) - ...
        innerRetinaRef.mosaic{1}.tonicDrive{cellNum};
    rLinSUVar(cellNum, 1, :) = rLinSU(cellNum, 1, :) - ...
        innerRetinaUn.mosaic{1}.tonicDrive{cellNum};

    rLinSUNew(cellNum, 1, :) = max(rLinVar(cellNum, 1, :)) * ...
        (rLinSUVar(cellNum, 1, :) / ...
        max(rLinSUVar(cellNum, 1, :))) + ...
        innerRetinaUn.mosaic{1}.tonicDrive{cellNum};
    % rLinSUNew{cellNum, 1, 1} = ...
    %    [rLinSUNew{cellNum, 1, 1} min(rLin{cellNum, 1})];

    % maxTemp = max(rLin{cellNum, 1});
    % maxSU = max(rLinSU{cellNum, 1});
    % rLinSUNew{cellNum, 1, 1} = rLinSU{cellNum, 1, 1} - (maxSU - maxTemp);
end

% figure;
% plot(rLinSUNew{cellNum, 1, 1}); hold on; plot(rLin{cellNum, 1, 1}, 'r');
innerRetinaUn.mosaic{1} = ...
    mosaicSet(innerRetinaUn.mosaic{1}, 'responseLinear', rLinSUNew);
innerRetinaUn.mosaic{1} = ...
    innerRetinaUn.mosaic{1}.mosaicSet('responseLinear', rLinSUNew);

end
