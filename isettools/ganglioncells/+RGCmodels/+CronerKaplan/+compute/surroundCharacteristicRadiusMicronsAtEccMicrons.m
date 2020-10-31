function radiusMicrons = surroundCharacteristicRadiusMicronsAtEccMicrons(eccMicrons, model)
    % method coould be the model fit by C&P to both M and P, or something
    % else

    eccDegs = RGCmodels.Watson.convert.rhoMMsToDegs(eccMicrons*1e-3);
    switch (model)
        case 'Croner_Kaplan_Fit_To_P_and_M_data'
            radiusDegs = RGCmodels.CronerKaplan.constants.surroundCharacteristicRadiusDegsAtEccentricity(eccDegs);
        otherwise
            error('Unknown model: ''%s''.', model);
    end
    radiusMicrons = RGCmodels.Watson.convert.sizeVisualDegsToSizeRetinalMicrons(radiusDegs);
end
