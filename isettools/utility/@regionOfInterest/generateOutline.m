function outline = generateOutline(obj)

    switch (obj.geometryStruct.shape)
        case 'rect'
            a = 0.5;
            b = 0.5;
            cosOutline = [-1 -1  1  1 -1]*obj.geometryStruct.width;
            sinOutline = [-1  1  1 -1 -1]*obj.geometryStruct.height;
        case 'ellipse'
            a = 0.5*obj.geometryStruct.minorAxisDiameter;
            b = 0.5*obj.geometryStruct.majorAxisDiameter;
            dTheta = 10;
            cosOutline = cosd(0:dTheta:360);
            sinOutline = sind(0:dTheta:360);
    end
    
    outline = struct(...
        'x', obj.geometryStruct.center(1) + a*cosOutline*cosd(obj.geometryStruct.rotation) - b*sinOutline*sind(obj.geometryStruct.rotation), ...
    	'y', obj.geometryStruct.center(2) + a*cosOutline*sind(obj.geometryStruct.rotation) + b*sinOutline*cosd(obj.geometryStruct.rotation) ...
        );
end