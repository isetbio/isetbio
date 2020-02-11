function theMeridianAngle = meridianAngle(obj, meridianName, meridianSpace, whichEye)
    
    % Validate input
    obj.validateMeridianName(meridianName);
    obj.validateMeridianSpace(meridianSpace);
    obj.validateEye(whichEye);

    switch (meridianName)
        case 'temporal meridian'
            if (strcmp(whichEye, 'right'))
                theMeridianAngle = 180;
            else
                theMeridianAngle = 0;
            end
        case 'superior meridian'
            theMeridianAngle = 90;
        case 'nasal meridian'
            if (strcmp(whichEye, 'right'))
                theMeridianAngle = 0;
            else
                theMeridianAngle = 180;
            end
        case 'inferior meridian'
            theMeridianAngle = 270;
        otherwise
            fprintf('Valid meridian names are: %s\n', keys(obj.meridianParams));
            error('Invalid meridian name: ''%s''.', meridian);
    end
    
    if (strcmp(meridianSpace, 'retinal'))
        theMeridianAngle = mod(theMeridianAngle + 180,360);
    end
end

