% Correct meridian color & legend based on which eye we are showing and whether we
% are labeling meridians in retinal or visual space
function [theLegend, theColor] = correctLegendAndColor(obj, theLegend, theColor, meridianName, displayRetinalMeridiansLegends, whichEye)
        
    if (displayRetinalMeridiansLegends)
        theLegend = retinalMeridianLegendFromRightEyeVisualFieldMeridianLegend(theLegend);
    end

    if ((strcmp(whichEye, 'left')) && (~displayRetinalMeridiansLegends))
        % Reverse the colors, since we swap the nasal/temporal angles in coneRFSpacingAndDensity()
        % to match the swapping done in coneDensityReadData for left vs right eye naso-temporal cone density data
        % This is so we can match the MeridianConventionsFigure, where
        % the right visual field is red-coded and the left is green coded
        if (contains(theLegend, 'nasal'))
            theColor = obj.meridianColors(strrep(meridianName,'nasal', 'temporal'));
        elseif (contains(theLegend, 'temporal'))
            theColor = obj.meridianColors(strrep(meridianName,'temporal', 'nasal'));
        end
    end
end

function theLegend = retinalMeridianLegendFromRightEyeVisualFieldMeridianLegend(theLegend)
     
    if (contains(theLegend, 'superior'))
        theLegend = strrep(theLegend, 'superior', 'inferior');
    elseif (contains(theLegend, 'inferior'))
        theLegend = strrep(theLegend, 'inferior', 'superior');
    elseif (contains(theLegend, 'nasal'))
        theLegend = strrep(theLegend, 'nasal', 'temporal');
    elseif (contains(theLegend, 'temporal'))
        theLegend = strrep(theLegend, 'temporal', 'nasal');
    end
    
    theLegend = strrep(theLegend, 'visual field', 'retina');
end