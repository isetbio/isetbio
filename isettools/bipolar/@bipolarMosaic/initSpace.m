function initSpace(obj, varargin)
% Initialize space required for the bipolar mosaic
%
% Syntax:
%
%   @bipolarMosaic.initSpace(varargin)
%
% Description:
%    Each bipolar mosaic takes its input from the cone mosaic, so that cell
%    locations are with respect to the spatial samples of the cone mosaic.
%    To compute the  spatial spread, we need to account for the cone
%    spacing. Thus, if the cones are spaced, say 2 um, and the bipolar RF
%    spans 5 samples, the spatial extent will be 2*5 um. 
%
% Input:
%
% Optional Key/Value Pairs:
%    eccentricity - placeholder for eccentricity, affects spread
%    spread       - placeholder for spread, in microns
%    spreadRatio  - placeholder for spreadRatio
%    stride       - placeholder for stride
%    ampCenter    - placeholder for ampCenter
%    ampSurround  - placeholder for ampSurround
%
% Notes:
%    * Parasol is synonymous with diffuse.
%    * TODO
%       - To compute the spread in microns from this specification, 
%         multiply the number of input samples by the spatial sample
%         separation of the cones in the mosaic (stored in the input slot).
%       - We will incorporate a function that changes the size of the
%         spread and support as a function of eccentricity. For now we just
%         put in some placeholder numbers. (Let's do better on this
%         explanation, BW).
%       - When the layer is deeper, however, we have to keep referring back
%         through multiple layers. This issue will be addressed in the
%         RGCLAYER, and then onward.
%       - We need to write simple utilities that convert from the spatial
%         units on the cone mosaic into spatial units on the retinal
%         surface (in um or mm). That will be first implemented in
%         bipolar.plot('mosaic'). But basically, to do this the units are
%         X*coneMosaic.patternSampleSize (we think). This doesn't deal with
%         the jittered cone mosaic yet, but kind of like this. (BW/JRG). 
% 
%
% References:
% * Field and Chichilnisky, 2010, Nature
%   <http://www.nature.com/nature/journal/v467/n7316/full/nature09424.html>
% * Dacey, Brainard, Lee, et al., Vision Research, 2000.
%   <http://www.cns.nyu.edu/~tony/vns/readings/dacey-etal-2000.pdf>
% * Size of the RF
% * Sampling density (stride) of the RF centers.
% * References and Built-in Bipolar Types
%       We have implemented five types of bipolar receptive fields, one
%       assigned to each of the big five RGC types. Each bipolar type has a
%       preferential cone selection rule. The critical rule is **no S-cones
%       for on/off parasol and on-midget**, as per the Chichilnisky primate
%       data (Field and Chichilnisky, 2010, Nature).
%   <http://www.nature.com/nature/journal/v467/n7316/full/nature09424.html>
% * N.B. Parasol is synonymous with diffuse.
% * The data for the support size is this passage from Dacey, Brainard,
%   Lee, et al., Vision Research, 2000, page 1808 bottom right.
% <http://www.cns.nyu.edu/~tony/vns/readings/dacey-etal-2000.pdf>
%   They write:
%         "The [spatial] frequency response was bandpass and well fit by a
%          difference of Gaussians receptive field model. (abstract)"
%
%         "For midget bipolar cells, it is known that at retinal
%          eccentricities up to 10 mm virtually all cells restrict
%          dendritic contact to single cones (Milam et al., 1993; Wassle et
%          al., 1994); this was confirmed for the cell whose light response
%          is illustrated in Fig. 4. B) Also see Boycott & Wassle, 1991,
%          (European Journal of Neuroscience), Table 1."
%   On page 1809:
%       Center/Surround gain ratio is about 1:1.3 (area under the curve)
%       Surround:Center diameter about 1:10 (Center:surround)
%       They seem to think that for ganglion cells the gain ratio is about
%       1:0.5 and the diameter ratio is between 1:2 and 1:5.
%
%       Likely the larger RF sizes measured physiological (Dacey et al.) vs
%       anatomically (B&W) reflect spread of signals among cones (via
%       direct gap junctions) and probably more important among cone
%       bipolars (via gap junctions with AII amacrine cells). - Fred
%

%% History:
% JRG/BW ISETBIO Team, 2015
%
%    10/19/17  jnm  Comments & Formatting

%% Parse inputs
p = inputParser;
p.KeepUnmatched = true;
p.addParameter('eccentricity', 0, @isscalar);
p.addParameter('spread', [], @isscalar);
p.addParameter('spreadRatio', 10, @isscalar);  % Dacey bipolar paper

p.addParameter('stride', [], @(x)(isempty(x) || isscalar(x)));
p.addParameter('ampCenter', 1, @(x)(isempty(x) || isscalar(x)));
p.addParameter('ampSurround', 1.3, @(x)(isempty(x) || isscalar(x)));

p.parse(varargin{:});

eccentricity = p.Results.eccentricity;
spread       = p.Results.spread;             % Spread center
stride       = p.Results.stride;
ampCenter    = p.Results.ampCenter;
ampSurround  = p.Results.ampSurround;
spreadRatio  = p.Results.spreadRatio;  % Surround spread / center spread

%% Select parameters for each cell type
% The spatial samples below (e.g. support and spread) are in units of
% samples on the cone mosaic. We can convert this to spatial units on the
% cone mosaic (microns) by multiplying by the cone spatial sampling
% distance. The cone mosaic is stored in the input slot of the bipolar
% mosaic.
switch obj.cellType
    
    case{'ondiffuse', 'offdiffuse'}
        %%%
        % Diffuse bipolars that carry parasol signals
        if isempty(spread)
            %%%
            % A functional rule, you could use this. 1 at the central
            % fovea and 3 cones at 40 deg. Linear.
            spread =  floor(1 + (2/40)*(eccentricity));
        end
        support = round(4*spread);    % Minimum spatial support
        %%%
        % Standard deviation of the Gaussian for the center, specified in
        % spatial samples on the input mosaic. Anywhere near the center
        % the input is basically 1 cone. Far in the periphery, it will be
        % seomthing else that we will have a function for, like the
        % support.
        %
        % We need an amplitude for these functions to be specified in the
        % object.
        obj.sRFcenter   = ampCenter*fspecial('gaussian', [support, ...
            support], spread);
        %%%
        % Calculate the spread of the surround
        spreadSurround = spread*spreadRatio;
        obj.sRFsurround = ampSurround*fspecial('gaussian', [support, ...
            support], spreadSurround);
            
    case {'onsbc'}
        %%%
        % Small bistratified cells - handle S-cone signals
        % Find reference from Fred/EJ.
        if isempty(spread)
            %%%
            % A functional rule, you could use this. 2 at the central
            % fovea and 5 cones at 40 deg. Linear.
            spread =  floor(2 + (3/40)*(eccentricity));
        end
        support = round(4*spread);    % Minimum spatial support
        %%%
        % Reference for very broad surround is in header
        obj.sRFcenter   = ampCenter*fspecial('gaussian', [support, ...
            support], spread);    % convolutional for now
        
        spreadSurround = spread*spreadRatio;
        obj.sRFsurround = ampSurround*fspecial('gaussian', [support, ...
            support], spreadSurround); % convolutional for now

    case{'onmidget', 'offmidget'}
        %%%
        % Midget bipolars to midget RGCs. Midgets restrict their centers
        % to 1 cone out to at least 10 mm in macaque, which is like 40 or
        % 50 deg.
        %
        % See the Psychtoolbox external routine RetinalMMToDegrees for help
        % with mm to degrees in different species.

        if isempty(spread)
            %%%
            % A functional rule, you could use this. 1 at the central
            % fovea and 3 cones at 40 deg. Linear.
            spread =  floor(1 + (2/40)*(eccentricity));
        end
        support = round(4*spread);    % Minimum spatial support

        % ecc = 0 mm yields 1x1 cone input to bp
        % ecc = 30 mm yields 3x3 cone input to bp
        % Support formula extrapolated from data in Dacey ... Lee, 1999 
        % @JRG to insert
        % support = max(minSupport, floor(1 + (2/10)*(eccentricity)));
        %%%
        % Standard deviation of the Gaussian for the center. Anywhere near
        % the center the input is basically 1 cone. Far in the periphery, 
        % it will be seomthing else that we will have a function for, like
        % the support.s
        obj.sRFcenter   = ampCenter*fspecial('gaussian', [support, ...
            support], spread); % convolutional for now
        %%%
        % Calculate the spread of the surround
        spreadSurround = spread*spreadRatio;
        obj.sRFsurround = ampSurround*fspecial('gaussian', [support, ...
            support], spreadSurround); % convolutional for now

end
%%%
% The bipolar RF center positions are stored with respect to the samples of
% the input layer (cone mosaic). The weights are also stored with respect
% to the input sample.
if isempty(stride), stride = round(spread); end
%%%
% Cone row and column positions, but centered around (0, 0).
% These should be spaced by an amount that is controlled by a parameter and
% reflects the size of the receptive field.
conemosaic = obj.input;
[X, Y] = meshgrid(1:stride:conemosaic.cols, 1:stride:conemosaic.rows);
X = X - mean(X(:)); Y = Y - mean(Y(:));
%%%
% Put them in the (row, col, X/Y) tensor.
obj.cellLocation = zeros(size(X, 1), size(X, 2), 2);
obj.cellLocation(:, :, 1) = X;
obj.cellLocation(:, :, 2) = Y;

end
