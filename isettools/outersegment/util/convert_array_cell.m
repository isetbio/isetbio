
function ConeCurrentSignalCell = convert_array_cell(params, sensor, os1);


cone_mosaic = sensorGet(sensor,'cone type');
[csz1, csz2] = size(cone_mosaic);
for cone_type = 2:4

    cone_locs = find(cone_mosaic==cone_type);
    ConeSignal_rs = reshape(os1.ConeCurrentSignal,[csz1*csz2],params.nsteps);
    ConeSignalFinal_rs(cone_locs,:) = ConeSignal_rs(cone_locs,:);
    ConeCurrentSignalCell{cone_type-1} = ConeSignal_rs(cone_locs,:);
        
end
