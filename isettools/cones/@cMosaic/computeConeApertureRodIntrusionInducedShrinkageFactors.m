function computeConeApertureRodIntrusionInducedShrinkageFactors(obj)

    % Compute radial eccDegs of all cones
    coneEccentricityDegs = sqrt((obj.coneRFpositionsDegs(:,1)).^2 + (obj.coneRFpositionsDegs(:,2)).^2);
   
    % Compute shrinkage at eccentricities
    obj.coneApertureRodIntrusionInducedShrinkageFactors = coneApertureShrinkageFactorsAtEccentricity(coneEccentricityDegs);
end

function shrinkageFactors = coneApertureShrinkageFactorsAtEccentricity(eccDegs)

    measurement = 'DirectFromCurcioPaper';
    %measurement = 'FromBanksPaper';

    if (strcmp(measurement, 'FromBanksPaper'))
        % Scanned data (eccentricity, osLength in microns) from Figure 1 (left panel)
        % Banks, Sekuler and Anderson (1991). Peripher spatial vision: limits
        % imposed by optics, photoreceptors and receptor pooling
        s = [ ...
            0.0000000000000  	0.57   0.57;
            2.042432662814E0	1.18   1.10;
            5.075542489714E0	1.81   1.33;
            9.920866641161E0	2.25   1.58;
            1.994845257032E1	3.02   1.63;
            4.016571693212E1	3.59   1.71;
        ];
    
      eccDegsRaw = s(:,1);
      coneSpacingRaw = s(:,2);
      innerSegmentDiameterRaw = s(:,3);

    else
      squaredMicrons = 38*26;
      sMMs = [ ...
          0    1 1;
          0.66 RGCmodels.Watson.convert.densityToSpacingForHexGrid(29/squaredMicrons) 5.9;
          1.35  RGCmodels.Watson.convert.densityToSpacingForHexGrid(18.5/squaredMicrons) 6.5;
          5.0  RGCmodels.Watson.convert.densityToSpacingForHexGrid(8/squaredMicrons) 7;
          8.0  RGCmodels.Watson.convert.densityToSpacingForHexGrid(6/squaredMicrons) 7.5;
          16.0 RGCmodels.Watson.convert.densityToSpacingForHexGrid(4/squaredMicrons) 8];
    
      eccMMsRaw = sMMs(:,1);
      coneSpacingRaw = sMMs(:,2);
      innerSegmentDiameterRaw = sMMs(:,3);
    
      eccDegsRaw = RGCmodels.Watson.convert.rhoMMsToDegs(eccMMsRaw);
    end


  interpolationMethod = 'pchip';
  shrinkageFactors = interp1(eccDegsRaw, innerSegmentDiameterRaw ./ coneSpacingRaw, eccDegs, interpolationMethod);
  

%   figure(10); hold on
%   plot(linspace(0,30,100), interp1(eccDegsRaw, innerSegmentDiameterRaw ./ coneSpacingRaw, linspace(0,30,100), 'pchip'), 'r-');
%   xlabel('ecc (degs)');
%   ylabel('aperture diameter/spacing')
%   set(gca, 'YLim', [0 1]);
%   pause
  shrinkageFactors = reshape(shrinkageFactors, [1 size(eccDegs,1)]);
end


