% Method to export cone positions, cone spacings, RGC RF positions
% and spacings, and an [Nx2] matrix  of all connection tuplets [mRGCIdx, coneIdx]
function  connectionsList = exportRFCenterWiring(obj, exportDataFileName)

    s = struct(...
        'mRGCsNum', size(obj.coneConnectivityMatrix,2), ...
        'conesNum', size(obj.coneConnectivityMatrix,1), ...
        'connections', []);

    
    idx = find(obj.coneConnectivityMatrix(:) > 0);
    connectionsNum = numel(idx);

    tic
    connectionsList = zeros(connectionsNum,2, 'uint32');
    connectionIndex = 0;
    noCenterInputs = 0;
    centerInputs = zeros(1,100);

    for RGCindex = 1:s.mRGCsNum
        connectivityVector = full(squeeze(obj.coneConnectivityMatrix(:, RGCindex)));
        inputConeIDs = find(connectivityVector > 0.01);
        if (numel(inputConeIDs) == 0)
            noCenterInputs = noCenterInputs+1;
        else
            centerInputs(numel(inputConeIDs)) = centerInputs(numel(inputConeIDs)) + 1;
        end
        for inputConeIndex = 1:numel(inputConeIDs)
            connectionIndex = connectionIndex + 1;
            connectionsList(connectionIndex,:) = uint32([RGCindex inputConeIDs(inputConeIndex)]);
        end
    end
    toc
    
    fprintf('There are %d RGC with 0 cones in their RF center.\n', noCenterInputs);
    for centerInputsNo = 1:100
        fprintf('There are %d RGC with %d cones in their RF center.\n', centerInputs(centerInputsNo), centerInputsNo);
    end
    
    if (~isempty(exportDataFileName))
        conePos = obj.inputConeMosaic.coneRFpositionsDegs;
        coneSpacings = obj.inputConeMosaic.coneRFspacingsDegs;
        coneTypes = obj.inputConeMosaic.coneTypes;
        RGCpos = obj.rgcRFpositionsDegs;
        RGCspacings = obj.rgcRFspacingsDegs;
        save(exportDataFileName, 'conePos', 'RGCpos', 'coneSpacings', 'coneTypes', 'RGCspacings', 'connectionsList');
    end
    
end
