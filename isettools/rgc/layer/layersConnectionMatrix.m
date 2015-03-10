function [connecMatrix, layer] = layersConnectionMatrix(layer)
% Calculate a lateral connection matrix between RGCs
%
%    [connecMatrix layer] = layersConnectionMatrix(layer)
%
% The connection strength falls off exponentially with distance from the
% RGC, plus (a) we add some noise, and we set a cutoff level below which
% connection strengths are set to 0.
%
% Here, we should give a physical meaning to the connection weight.
% 
% The matrix connecMatrix is a cell array with the following structure:
%
%   connecMatrix(ii,1) = [i1 i2 i3 ... in] corresponding to the indices of 
% 
% the cells affected by the ii-th cell
% 
%   connecMatrix(ii,2) = [w1 w2 w3 ... wn]: corresponding connection weights
% 
% We need an explanation of what the weights represent for each connection.
% 
% This function is more memory friendly than rgcLatCreateFull which
% creates a 2*nRGC*nRGC connection matrix.
%
% The weights depend on the distance function (f): f(wRange,d) + nRange*rand(),
% where f is a determistic function specified in param.distanceFunction, it
% is a function of wRange and the distance, d
%
% The function is assuming there is no link between the layers.
%
% TODO: 
% 1. The logic of building these connections needs to be explained
% better in the code below. 
% 2. There are a lot of loops in the function, make it faster 
% (low usefulness as done only once => low priority)
%
%
%Examples:
% First, create the parameters
%   rgcP = rgcParameters;
%   rgcP.addOnLayer; rgcP.set('data',rand(128,128,10));
%   rgcP = rgcComputeLayers(rgcP);
%
%   rgcP.wRange = 0; rgcP.nRange = 0;   % No lateral connections
%   [latMatrix rgcP] = rgcCreateConnectionMatrix(rgcP);  
%   max(latMatrix(:,:,2))
%
%   rgcP.wRange = 100; rgcP.nRange = 0; % Massive, no noise connections
%   [latMatrix rgcP] = rgcCreateConnectionMatrix(rgcP);  
%   figure(1); hist(latMatrix,50);
%   max(latMatrix(:,:,2)), min(latMatrix(:,:,2))
%
%   rgcP.wRange = 100; rgcP.nRange = 10; % Massive, noisy connections
%   [latMatrix rgcP] = rgcCreateConnectionMatrix(rgcP);  
%   figure(2); hist(latMatrix,50);
%   max(latMatrix(:,:,2)), min(latMatrix(:,:,2))
%
% (c) Synapse project Stanford 2009

%% Check network input
if (notDefined('layer') || ~isequal(class(layer),'rgcLayer'))
    error('There is no argument layer belonging to the rgcLayer class');
end;

%% Check if needed to be computed:
gridSize = layer.get('gridSize');
nRGC     = gridSize(1)*gridSize(2);

% currentConnec is a cell array with two dimensions.
%
% The first dimension represents which RGC.
% For each RGC we have two entries.  The first defines which (random
% number) of cells are connected to this one.  Suppose there are C
% connected. 
% The second cell entry defines the weights of the connections between this
% cell and others???  I am not sure why there are always (C,100) values in
% the 2nd slot.
%
currentConnec = layer.get('currentConnec');
if ~isempty(currentConnec)
    connecMatrix = currentConnec; % We stop here
else
    
    %% Extract connection parameters 
    connecMatrix = cell(nRGC,2);

    % wRange is the parameter of the distance function
    wRange  = layer.get('wRange');

    % nRange is the noise variance for the connections
    nRange  = layer.get('nRange');

    % The threshold for setting weights to 0
    cutoff = layer.get('cutoff'); 

    % Inline function of distance to weight
    distanceFunction = layer.parent.get('distanceFunction');

    %% Compute the maximal size of the cube influenced by a cell 
    %  In general there is no link between layers, so the cube is a square
    %
    % We need to have (nRange << cutoff) otherwise the matrix will not be
    % sparse. So we impose (nRange <= cutoff/2)
    if (nRange > cutoff/2)
        error('For memory reasons, we impose nRange <= cutoff/2, To use your settings use rgcLatCreateFull');
    end

    threshold = cutoff-nRange;
    % if distanceFunction(d,wRange)<threshold then there can not be a
    % connection between the 2 cells.

    % Grid distance between 2 cells
    cellSpacing = layer.get('cellSpacing');
    cubeSide=0;

    while distanceFunction(cubeSide*cellSpacing,wRange) >= threshold
        cubeSide = cubeSide+1;
    end
    cubeSide = max(cubeSide-1,0);

    if cubeSide>3
        warning('Settings of the distance function and cellSpacing are such that each cell has an impact to a distance bigger than 3. The coupling step is likely to take a really long time.');
        yes = input('Do you want to stop? ');
        if yes 
            error('distance settings - user stopped');
        end
    end

    if cubeSide > 0

        %% Calculate the connection between each neuron in the cube and the center
        middle = cubeSide+1;
        theCube = zeros(2*cubeSide+1);
        for ii = 1:middle
            for jj = 1:ii
                di = middle-ii;
                dj = middle-jj;
                d = distanceFunction(cellSpacing*sqrt(di*di+dj*dj),wRange);
                i1 = ii;
                j1 = jj;
                theCube(i1,j1) = d;
                theCube(2*cubeSide+1-i1+1,j1) = d;
                theCube(i1,2*cubeSide+1-j1+1) = d;
                theCube(2*cubeSide+1-i1+1,2*cubeSide+1-j1+1) = d;
                i2 = jj;
                j2 = ii;
                theCube(i2,j2) = d;
                theCube(2*cubeSide+1-i2+1,j2) = d;
                theCube(i2,2*cubeSide+1-j2+1) = d;
                theCube(2*cubeSide+1-i2+1,2*cubeSide+1-j2+1) = d;
            end
        end

        %% Build connection matrix
        C = zeros(nRGC,(2*cubeSide+1)*(2*cubeSide+1),2);

        % Looping through all the cells to build C
        for ii = 1:gridSize(1);
            for jj = 1:gridSize(2)
                indice = (ii-1)*gridSize(2)+jj;
                n = 0;
                for aa = -cubeSide:cubeSide
                    for bb = -cubeSide:cubeSide
                        % Not treating self connections
                        if (aa || bb)
                            iii = ii + aa;
                            jjj = jj + bb;
                            indice2 = (iii-1)*gridSize(2)+jjj;
                            if (iii>0 && iii<=gridSize(1) && jjj>0 && jjj<=gridSize(2))
                               n = n+1; 
                               C(indice,n,1) = indice2;
                               C(indice,n,2) = theCube(cubeSide+1+aa,cubeSide+1+bb);
                            end
                        end
                    end
                end
            end
        end

        % Adding the noise to the connection weights
        C(:,:,2) = C(:,:,2)+ nRange*rand(nRGC,(2*cubeSide+1)*(2*cubeSide+1),1);

        %% Finish connection matrix

        C2 = squeeze(C(:,:,2));
        % Apply the cutoff value
        C2(C2 < cutoff) = 0;

        C(:,:,2) = C2;

        %assuming they are all identical:
        cpTR  = layer.get('cpTR');   % cp probably means coupling; TR??
        for ii = 1:size(C,1); % for each rgc
            b = squeeze(C(ii,:,:));
            a = sortrows(b,-2);
            cpW = nonzeros(a(:,2)); % non-zero connexion values, sorted, biggest first
            jj = a(:,1); 
            connecMatrix{ii,1} = jj(1:length(cpW)); % indices of connected neurons (non-zero connexions)
            connecMatrix{ii,2} = cpW*cpTR;          % connexion values, multiplied by cpTR
        end
    end

    layer.set('currentConnec',connecMatrix); 

end
