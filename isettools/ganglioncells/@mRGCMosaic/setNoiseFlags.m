function setNoiseFlags(obj, varargin)
    % Parse input
    p = inputParser;
    p.addParameter('noiseFlag', [], @(x)(isempty(x))||(ismember(x, {'none', 'frozen', 'random'})));
    p.addParameter('randomSeed', [], @(x)(isempty(x))||(isinteger(x)&&(numel(x)==1)));
    p.addParameter('inputConeMosaicNoiseFlag', [], @(x)(isempty(x))||(ismember(x, {'none', 'frozen', 'random'})));
    p.addParameter('inputConeMosaicRandomSeed', [], @(x)(isempty(x))||(isinteger(x)&&(numel(x)==1)));
    p.parse(varargin{:});
    

    if (~isempty(p.Results.noiseFlag))
        switch (p.Results.noiseFlag)
            case 'none'
                % Do nothing with the randomSeed
    
            case 'random'
                % Make a noise seed for the input cone mosaic and store it.  One in a million.
                obj.randomSeed = randi(1e6);
                rng(obj.randomSeed);
    
            case 'frozen'
                if (~isempty(p.Results.randomSeed))
                    % Check that a random seed  was provided and set the rng
                    obj.randomSeed = p.Results.randomSeed;
                    rng(obj.randomSeed);
                else
                    error('Frozen noise but no seed provided.');
                end
        end
        obj.noiseFlag = p.Results.noiseFlag;
    end

    if (~isempty(p.Results.inputConeMosaicNoiseFlag))
        switch (p.Results.inputConeMosaicNoiseFlag)
            case 'none'
                % Do nothing with the input cone mosaic randomSeed
            case 'random'
                % Make a noise seed for the input cone mosaic and store it.  One in a million.
                obj.inputConeMosaic.randomSeed = randi(1e6);
    
            case 'frozen'
                if (~isempty(p.Results.inputConeMosaicRandomSeed))
                    % Check that a random seed for the input cone mosaic was provided
                    obj.inputConeMosaic.randomSeed = p.Results.inputConeMosaicRandomSeed;
                else
                    error('Frozen noise but no seed provided for the input cone mosaic.');
                end
        end
        obj.inputConeMosaic.noiseFlag = p.Results.inputConeMosaicNoiseFlag;
    end

end
