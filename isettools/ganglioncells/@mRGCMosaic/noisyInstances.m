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
        case 'frozen'
            rng(seed);
        case 'random'
            rng('shuffle');
    end

    % vMembrane additive noise
    fprintf('Computing noisy mRGC response instances with vMembrane noise std = %2.3f\n', obj.vMembraneGaussianNoiseSigma);
    noisyMRGCresponseInstances = noiseFreeMRGCresponses + ...
                obj.vMembraneGaussianNoiseSigma * randn(size(noiseFreeMRGCresponses));

end