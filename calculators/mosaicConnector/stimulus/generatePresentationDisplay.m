function presentationDisplay = generatePresentationDisplay(varargin)
    % Parse input
    p = inputParser;
    p.addParameter('viewingDistance', 0.57, @isscalar);  % 1 cm at 57 cm occupies 1 degree of visual space
    p.addParameter('bitDepth', 25, @isscalar);
    p.addParameter('gammaTable', 'linear', @(x)(ismember(x,{'linear', 'nonlinear'})));
    p.addParameter('wave', 400:10:750, @isnumeric);
    p.parse(varargin{:});
            
    % Generate presentation display
    presentationDisplay = displayCreate('LCD-Apple', ...
        'wave', p.Results.wave, ...
        'viewing distance',p.Results.viewingDistance);
    bitDepth = p.Results.bitDepth;
    gammaTable = p.Results.gammaTable;
    if (strcmp(gammaTable, 'linear'))
        N = 2^bitDepth;
        gTable = repmat(linspace(0, 1, N), 3, 1)';
    else
        x = linspace(0,1,2^bitDepth );
        gTable = x(:).^2.1;
        gTable = repmat(gTable, [1,3]);
    end
    presentationDisplay = displaySet(presentationDisplay, 'gTable', gTable);
end
