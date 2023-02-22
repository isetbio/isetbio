function  opticsPositionDegs = selectOpticsPosition(theMidgetRGCmosaic)

    opticsPositionDegs = [];
    while (numel(opticsPositionDegs) ~= 2)
        opticsPositionDegs = input('\nEnter the optics position to use for computing retinal stimulus images ([x y]): ');
    end
end
