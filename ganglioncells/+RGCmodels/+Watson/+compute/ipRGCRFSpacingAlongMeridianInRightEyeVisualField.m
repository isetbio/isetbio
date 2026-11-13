function [ipRGCSpacingDegs, ipRGCSpacingMMs, ipRGCDensityDegs2, ipRGCDensityMMs2, ...
    scannedData] = ipRGCRFSpacingAlongMeridianInRightEyeVisualField(eccDegs, rightEyeVisualFieldMeridianName)
% Compute ipRGC receptive field spacing&density for desired positions (specified in visual degrees) along a meridian on the right eye visual field
% Note: This is for the combination of all ipRGCs subtypes

     % Ensure all eccentricities are positive
    assert(all(eccDegs>=0), sprintf('Eccentrity vector must not include any negative values'));


    % Convert eccDegs to eccMMs so that we don't go through the baked-in
    % ecc-independent 0.3 mm / deg conversion in coneDensityReadData()
    eccMM = RGCmodels.Watson.convert.rhoDegsToMMs(eccDegs);

    [scannedData.eccMMs, scannedData.densityCellsPerMMs2] = scannedDataFromDacey2005Figure1G([]);

    [~,ipRGCDensityMMs2] = scannedDataFromDacey2005Figure1G(eccMM);
    

    % ipRGC density to spacing in MMs
    ipRGCSpacingMMs = RGCmodels.Watson.convert.densityToSpacingForHexGrid(ipRGCDensityMMs2);
    
    % Convert density per mm^2 to density per deg^2
    ipRGCDensityDegs2 = RGCmodels.Watson.convert.densityMMs2ToDegs2(ipRGCSpacingMMs, eccDegs);
    
    % Convert density to spacing in degs
    ipRGCSpacingDegs = RGCmodels.Watson.convert.densityToSpacingForHexGrid(ipRGCDensityDegs2);
end

function [eccMMs, densityCellsPerMMs2] = scannedDataFromDacey2005Figure1G(desiredEccMMs)
d = [...
    2.041884816754E0	2.106205250597E1; ...
    3.047120418848E0	1.545346062053E1; ...
    4.052356020942E0	1.199284009547E1; ...
    5.026178010471E0	8.233890214797E0; ...
    6.031413612565E0	5.966587112172E0; ...
    7.036649214660E0	8.412887828162E0; ...
    8.041884816754E0	5.787589498807E0; ...
    9.078534031414E0	6.026252983294E0; ...
    1.005235602094E1	5.071599045346E0; ...
    1.099476439791E1	5.906921241050E0; ...
    1.196858638743E1	5.011933174224E0; ...
    1.300523560209E1	4.773269689737E0; ...
    1.401047120419E1	7.458233890215E0; ...
    1.507853403141E1	4.713603818616E0; ...
    1.598952879581E1	6.264916467780E0; ...
    1.696335078534E1	4.057279236277E0];

    eccMMs = d(:,1);
    densityCellsPerMMs2 = d(:,2);


    if (~isempty(desiredEccMMs))
        % Fit an exponential to the raw data
        [eccMMs, densityCellsPerMMs2] = fitExponentialWithPlateau(eccMMs, densityCellsPerMMs2, desiredEccMMs);
    
        densityCellsPerMMs2 = reshape(densityCellsPerMMs2, size(desiredEccMMs));
    end
end


function [xPred, yPred] = fitExponentialWithPlateau(x,y, xf)
    x = x(:);
    y = y(:);
    xf = xf(:);
    %% Model in terms of q = [A log10(tau) C]
    model = @(q, x) q(1)*exp(-x ./ 10.^q(2)) + q(3);

    %% Multi-start over tau; A and C are started from a linear fit for each tau
    xr      = max(x) - min(x);
    tauGrid = logspace(log10(0.05*xr), log10(2*xr), 12);     % do not start beyond 2x the x range
    best    = Inf;
    for k = 1:numel(tauGrid)
        ac = [exp(-x/tauGrid(k)), ones(size(x))] \ y;        % [A; C] for this tau
        q0 = [ac(1), log10(tauGrid(k)), ac(2)];
        try
            m = fitnlm(x, y, model, q0);
            if m.SSE < best, best = m.SSE; mdl = m; end
        catch
            % skip starts that fail to converge
        end
    end

    %% Results (tau reported in linear units)
    q   = mdl.Coefficients.Estimate';
    A   = q(1);  tau = 10^q(2);  C = q(3);
    ci  = coefCI(mdl);                                       % CIs in q-space
    tauCI = 10.^ci(2,:);                                     % monotone map, so endpoints transform directly
    
    if (1==2)
        fprintf('A   = %.4f  [%.4f, %.4f]\n', A,  ci(1,1), ci(1,2));
        fprintf('tau = %.4f  [%.4f, %.4f]\n', tau, tauCI(1), tauCI(2));
        fprintf('C   = %.4f  [%.4f, %.4f]\n', C,  ci(3,1), ci(3,2));
        fprintf('R^2 = %.5f,  RMSE = %.4f\n', mdl.Rsquared.Ordinary, mdl.RMSE);
    end
    
    if tau > 3*xr
        warning('tau (%.3g) is much larger than the x range (%.3g): the fit is nearly linear.', tau, xr);
    end

    yf = predict(mdl, xf);
    xPred = xf;
    yPred = yf;
end



function fromMRGC
    % Compute total RGC receptive field spacing&density
    [~, ~, totalrgcDensityDegs2, totalrgcDensityMMs2] = ...
        RGCmodels.Watson.compute.totalRGCRFSpacingAlongMeridianInRightEyeVisualField(eccDegs, rightEyeVisualFieldMeridianName);
    
    % Multiply with variation of midget RGC/totalRGC ratio with eccentricity
    rgcDensityDegs2 = RGCmodels.Watson.constants.midgetRGCFractionEccVariation(eccDegs) .* totalrgcDensityDegs2;
    rgcDensityMMs2  = RGCmodels.Watson.constants.midgetRGCFractionEccVariation(eccDegs) .* totalrgcDensityMMs2;
    
    % Convert density to spacing in degs
    rgcSpacingDegs = RGCmodels.Watson.convert.densityToSpacingForHexGrid(rgcDensityDegs2); 
    rgcSpacingMMs = RGCmodels.Watson.convert.densityToSpacingForHexGrid(rgcDensityMMs2); 
end


