function noisyMRGCresponseInstances = noisyInstances(obj, noiseFreeMRGCresponses, varargin)
    p = inputParser;
    p.addRequired('noiseFreeMRGCresponses', @isnumeric);
    p.addParameter('seed', 1, @isnumeric);
    p.parse(noiseFreeMRGCresponses, varargin{:});

    seed = p.Results.seed;

    if (strcmp(obj.noiseFlag, 'none'))
        noisyMRGCresponseInstances = [];
        return;
    end

    % Set up RNG depending on noiseFlag
    switch obj.noiseFlag
        case 'none'
            noiseGain = 0.0;
        case 'frozen'
            rng(seed);
            noiseGain = 1.0;
        case 'random'
            rng('shuffle');
            noiseGain = 1.0;
    end

    % vMembrane additive noise
    fprintf('Computing noisy mRGC response instances with vMembrane noise std = %2.3f\n', noiseGain * obj.vMembraneGaussianNoiseSigma);
    noisyMRGCresponseInstances = noiseFreeMRGCresponses + ...
                noiseGain * obj.vMembraneGaussianNoiseSigma * randn(size(noiseFreeMRGCresponses));

end