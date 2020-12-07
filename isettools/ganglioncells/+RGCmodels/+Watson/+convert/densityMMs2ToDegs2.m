function [densityPerDeg2, alpha] = densityMMs2ToDegs2(densityPerMM2, eccDegs)
% Convert retinal area mm^2 -> deg^2 as a function of 
% eccentricity in degs (Equation A7)

   % Compute mmSquaredPerDegSquared conversion factor alpha (Equation A7)
   % for the passed eccentricities (specified in degs)
   absEccDegs = abs(eccDegs);
   alpha = 0.0752...
         + 5.846 * 1e-5 * absEccDegs ...
         - 1.064 * 1e-5 * absEccDegs.^2 ...
         + 4.116 * 1e-8 * absEccDegs.^3;
        
    densityPerDeg2 = alpha .* densityPerMM2;       
end