function [isetbioAngle, whichEye, rightEyeRetinalMeridianName] = isetbioRetinalAngleForWatsonMeridian(obj,rightEyeVisualFieldMeridianName)
% Return ISETBio retinal angle for Watson's meridians (specified in visual field of the right eye)

    % Validate the meridianName
    obj.validateMeridianName(rightEyeVisualFieldMeridianName);
    
    % Return isetbio angle in the right eye
    whichEye = 'right';
    
    switch (rightEyeVisualFieldMeridianName)
        case obj.enumeratedMeridianNames{1}
            % Watson's temporal visual meridian in RE is ISETBio's nasal retina in the RE (180 degs) 
            isetbioAngle = 180;
            rightEyeRetinalMeridianName = 'nasal retina';
            
        case obj.enumeratedMeridianNames{2}
            % Watson's SUPERIOR visual meridian is ISETBio's INFERIOR retina (270 degs) 
            isetbioAngle = 90;
            rightEyeRetinalMeridianName = 'inferior retina';
            
        case obj.enumeratedMeridianNames{3}
            % Watson's NASAL visual meridian in RE is ISETBio's TEMPORAL retina in the RE (0 degs) 
            isetbioAngle = 0;
            rightEyeRetinalMeridianName = 'temporal retina';
            
        case obj.enumeratedMeridianNames{4}
            % Watson's INFERIOR visual meridian is ISETBio's SUPERIOR retina in the LE (90 degs) 
            isetbioAngle = 270;
            rightEyeRetinalMeridianName = 'superior retina';
    end
            
end
