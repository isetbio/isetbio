function noisyMRGCresponseInstances = noisyInstances(obj, noiseFreeMRGCresponses, varargin)
    p = inputParser;
    p.addRequired('noiseFreeMRGCresponses', @isnumeric);
    p.addParameter('noiseFlag', [], @(x)(isempty(x) | ismember(x, {'random', 'frozen','donotset','none'})));
    p.addParameter('seed', 1, @isnumeric);
    p.parse(noiseFreeMRGCresponses, varargin{:});

    seed = p.Results.seed;

    if (strcmp(obj.noiseFlag, 'none'))
        noisyMRGCresponseInstances = [];
        return;
    end

    % Set up RNG depending on noiseFlag.  This originally respected the
    % value in obj.noiseFlag.  We now pass a key/value pair that can
    % override use of that value, but were nervous about breaking existing
    % code simply by always using what is passed here.  So we default the
    % new key/value pair to empty and respect obj.noiseFlag if
    % p.Results.noiseFlag is empty.
    %
    % We noticed that the cMosaic noisyInstances method accepts a key/value
    % pair and that its value totally overrides its obj.noiseFlag, which is
    % either an error introduced inadvertently at some point, or strange
    % design.
    %
    % We can't currently set obj.noiseFlag to 'donotset' on object
    % creation, but we've left the case for it in case we ever do.
    %
    % At some moment when we have a little time, we should rationalize all
    % of this and deal with any backwards incompatibilities we introduce at
    % that point.
    if (isempth(p.Results.noiseFlag))
        switch obj.noiseFlag
            case 'none'
                noiseGain = 0.0;
            case 'frozen'
                rng(seed);
                noiseGain = 1.0;
            case 'random'
                rng('shuffle');
                noiseGain = 1.0;
            case 'donotset'
                noiseGain = 1.0;
        end
    else
        switch p.Results.noiseFlag
            case 'none'
                noiseGain = 0.0;
            case 'frozen'
                rng(seed);
                noiseGain = 1.0;
            case 'random'
                rng('shuffle');
                noiseGain = 1.0;
            case 'donotset'
                noiseGain = 1.0;
        end
    end

    % vMembrane additive noise
    fprintf('Computing noisy mRGC response instances with vMembrane noise std = %2.3f\n', noiseGain * obj.vMembraneGaussianNoiseSigma);
    noisyMRGCresponseInstances = noiseFreeMRGCresponses + ...
                noiseGain * obj.vMembraneGaussianNoiseSigma * randn(size(noiseFreeMRGCresponses));

end