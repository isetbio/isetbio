function outline = generateOutline(obj)

    switch (obj.geometryStruct.shape)
        case 'rect'
            a = 0.5;
            b = 0.5;
            cosOutline = [-1 -1  1  1 -1]*obj.geometryStruct.width;
            sinOutline = [-1  1  1 -1 -1]*obj.geometryStruct.height;

            outline = struct(...
                'x', obj.geometryStruct.center(1) + a*cosOutline*cosd(obj.geometryStruct.rotation) - b*sinOutline*sind(obj.geometryStruct.rotation), ...
                'y', obj.geometryStruct.center(2) + a*cosOutline*sind(obj.geometryStruct.rotation) + b*sinOutline*cosd(obj.geometryStruct.rotation) ...
                );
    
        case 'ellipse'
            a = 0.5*obj.geometryStruct.minorAxisDiameter;
            b = 0.5*obj.geometryStruct.majorAxisDiameter;
            dTheta = 10;
            cosOutline = cosd(0:dTheta:360);
            sinOutline = sind(0:dTheta:360);
            
            outline = struct(...
                'x', obj.geometryStruct.center(1) + a*cosOutline*cosd(obj.geometryStruct.rotation) - b*sinOutline*sind(obj.geometryStruct.rotation), ...
                'y', obj.geometryStruct.center(2) + a*cosOutline*sind(obj.geometryStruct.rotation) + b*sinOutline*cosd(obj.geometryStruct.rotation) ...
                );
            
        case 'line'
            % obj.geometryStruct
            rotation = atan2d((obj.geometryStruct.to(2)-obj.geometryStruct.from(2)),...
                              (obj.geometryStruct.to(1)-obj.geometryStruct.from(1)));
            width = obj.geometryStruct.thickness;
            length = sqrt(sum((obj.geometryStruct.from - obj.geometryStruct.to).^2));
            center = (obj.geometryStruct.from + obj.geometryStruct.to)/2;
            a = 0.5;
            b = 0.5;
            cosOutline = [-1 -1  1  1 -1]*length;
            sinOutline = [-1  1  1 -1 -1]*width;
            
            outline = struct(...
                'x', center(1) + a*cosOutline*cosd(rotation) - b*sinOutline*sind(rotation), ...
                'y', center(2) + a*cosOutline*sind(rotation) + b*sinOutline*cosd(rotation) ...
                );
    end
    
    
end