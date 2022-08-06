function  [semiAxes, rfCenter] = fillRFoutline(theAxes, C, zLevels, whichLevelsToContour, fitEllipse, faceAlpha, edgeAlpha )
    k = 1;
    contoursNum = 0;
    while k < size(C,2)
        level = C(1,k);
        points = C(2,k);
        if  (level ~= zLevels(whichLevelsToContour(1))) && (~ismember(level, zLevels(whichLevelsToContour)))
            % skip this contour
            k = k+points+1;
            continue;
        end
        
        xRGCEnsembleOutline = C(1,k+(1:points));
        yRGCEnsembleOutline = C(2,k+(1:points));
       
        if (fitEllipse)
            [xRGCEnsembleOutline,  yRGCEnsembleOutline, ...
                semiAxes, rfCenter, noFit] = fitEllipseToContour(xRGCEnsembleOutline,  yRGCEnsembleOutline);
        else
            semiAxes = [nan nan];
            rfCenter = [nan nan];
        end
        
        faceColor = [0.6 0.6 0.5]-level*0.05;
        edgeColor = [0.2 0.2 0.2];
        patchContour(theAxes, xRGCEnsembleOutline, yRGCEnsembleOutline, faceColor, edgeColor, faceAlpha, edgeAlpha);

        k = k+points+1;
        contoursNum = contoursNum + 1;
    end
end

function patchContour(theAxes, xRGCEnsembleOutline, yRGCEnsembleOutline, faceColor, edgeColor, faceAlpha, edgeAlpha)
    v = [xRGCEnsembleOutline(:) yRGCEnsembleOutline(:)];
    f = 1:numel(xRGCEnsembleOutline);
    patch(theAxes, 'Faces', f, 'Vertices', v, 'FaceColor', faceColor, ...
            'FaceAlpha', faceAlpha, 'EdgeColor', edgeColor, ... 
           'EdgeAlpha', edgeAlpha, 'LineWidth', 1.5);
end
