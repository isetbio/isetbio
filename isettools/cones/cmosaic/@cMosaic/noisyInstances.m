function noisyAbsorptionInstances = noisyInstances(meanAbsorptions, varargin)
    
    p = inputParser;
    p.addRequired('meanAbsorptions', @isnumeric);
    p.addParameter('noiseFlag', 'random', @(x)(ismember(x, {'random', 'frozen'})));
    p.addParameter('seed', -1, @isnumeric);
    p.parse(meanAbsorptions, varargin{:});

    seed = p.Results.seed;
    noiseFlag = p.Results.noiseFlag;

    % Set up RNG depending on noiseFlag
    %
    % rng('shuffle') is very slow.  Avoid unless
    % we really need to call it.
    if (seed ~= -1)
        rng(seed);
    else
        switch noiseFlag
            case 'frozen'
                if (seed == -1)
                    rng(1);
                else
                    rng(seed);
                end
            case 'random'
                rng('shuffle');
        end
    end

    % Gaussian approximation
    noise = sqrt(meanAbsorptions) .* randn(size(meanAbsorptions));
    noisyAbsorptionInstances = round(meanAbsorptions + noise);
    
    % For small counts (< 25) create a true Poisson sample.
    poissonCriterion = 25;
    idx = find(meanAbsorptions(:) < poissonCriterion);
    if ~isempty(idx)
        v = meanAbsorptions(meanAbsorptions < poissonCriterion);
        % Poisson samples
        noisyAbsorptionInstances(idx) = iePoisson(v, 'noiseFlag', noiseFlag, 'seed', seed);
    end

    % Can not have less than 0 absorptions
    noisyAbsorptionInstances(noisyAbsorptionInstances< 0) = 0;
end

