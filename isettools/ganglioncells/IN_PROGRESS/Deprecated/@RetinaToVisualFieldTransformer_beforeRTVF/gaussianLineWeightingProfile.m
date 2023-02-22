function theLineWeightingProfile = gaussianLineWeightingProfile(params, spatialSupport)
    % Retrieve params
    gain = params(1);
    xo = params(2);
    Rc = params(3);
    theLineWeightingProfile = gain * Rc * sqrt(pi) * exp(-((spatialSupport-xo)/Rc).^2);
end