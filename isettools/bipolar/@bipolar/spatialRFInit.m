function spatialRFInit(obj,varargin)
% Build spatial receptive field
%
% The data for the size of the support is based off of this passage
% from Dacey, Brainard, Lee, et al., Vision Research, 2000, page
% 1808 bottom right.
% (http://www.cns.nyu.edu/~tony/vns/readings/dacey-etal-2000.pdf)
% They write:
%
%  "The frequency response was bandpass and well fit by a difference of
%  Gaussians receptive field model. (abstract)"
%
%  "For midget bipolar cells, it is known that at retinal eccentricities up
%  to 10 mm virtually all cells restrict dendritic contact to single cones
%  (Milam et al., 1993; Wassle et al., 1994); this was confirmed for the
%  cell whose light response is illustrated in Fig. 4. B) Also see Boycott
%  & Wassle, 1991,  (European Journal of Neuroscience), Table 1."
%
%  On page 1809:
%   Center/Surround gain ratio is about 1:1.3 (area under the curve)
%   Surround:Center diameter about 1:10 (Center:surround)
%   They seem to think that for ganglion cells the gain ratio is about
%   1:0.5 and the diameter ratio is between 1:2 and 1:5.
%
% Likely the larger RF sizes measured physiological (Dacey et al.)
% vs anatomically (B&W) reflect spread of signals among cones (via
% direct gap junctions) and probably more important among cone
% bipolars (via gap junctions with AII amacrine cells). - Fred
%
% @JRG included  preferential cone selection to bipolar rfs. These
% set the basic parameters of the spatial receptive fields. There
% are some detailed modifications about the type of cone inputs (no
% S-cones to the on midget or the diffuse systems).
%
% We will also incorporate a function that changes the size of the
% spread and support as a function of eccentricity.  For now we
% just put in some placeholder numbers.

% These numbers don't make sense to BW at this time.  We
% need to write a script showing how big they are with
% respect to the cone mosaic, and we need to check how they
% vary with eccentricity.  Comparing with the curves in the
% cited data would be best.

%%
p = inputParser;
p.addParameter('eccentricity',0,@isscalar);

% For the future.  We don't have multiple mosaics yet.

p.parse(varargin{:});

eccentricity = p.Results.eccentricity;

%%
switch obj.cellType
    
    case{'ondiffuse','offdiffuse','onparasol','offparasol'}
        % Diffuse bipolars that carry parasol signals
        % ecc = 0 mm yields 2x2 cone input to bp
        % ecc = 30 mm yields 5x5 cone input to bp
        
        minSupport = 12;   % Minimum spatial support

        % BW, screwing around.  Just made spatial spread up here.
        % Support formula extrapolated from data in Dacey ... Lee, 1999 @JRG to insert
        support = max(minSupport,floor(2 + (3/10)*(eccentricity)));
        
        % Standard deviation of the Gaussian for the center.  Anywhere near
        % the center the input is basically 1 cone.  Far in the periphery,
        % it will be seomthing else that we will have a function for, like
        % the support.s
        spread = 1;  
        
        obj.sRFcenter   = fspecial('gaussian',[support, support],spread);  
        obj.sRFsurround = 1.3*fspecial('gaussian',[support,support], 1.3*spread);
        
        % How do we denote the spatial sample positions?  We refer it to
        % the positions of the cones in the coneMosaic.  The inputs to the
        % bipolars are actual cones, so we think it's OK to define the RF
        % in terms of weights from actual cones into the bipolar.
        %
        % As a utility, we want to be able to show this in spatial units on
        % the retinal surface (in um or mm).  That will be first
        % implemented in bipolar.plot('mosaic').
        obj.cellLocation = [];
        
    case {'onsbc'}
        minSupport = 15;    % Minimum spatial support

        % Small bistratified cells - handle S-cone signals
        
        % Needs to be checked and thought through some more @JRG
        % for this particular cell type.
        support = max(minSupport,floor(2 + (3/10)*(eccentricity)));
        
        spread = 3;  % Standard deviation of the Gaussian - will be a function
        rfCenterBig   = fspecial('gaussian',[support,support],spread); % convolutional for now
        rfSurroundBig = fspecial('gaussian',[support,support],10*spread); % convolutional for now
        
        obj.sRFcenter = rfCenterBig(:,:);
        obj.sRFsurround = rfSurroundBig(:,:);
        
        % How do spatially sample the positions?
        obj.cellLocation= [];
        
    case{'onmidget','offmidget'}
        % Midget bipolars to midget RGCs
        
        minSupport = 7;    % Minimum spatial support

        % ecc = 0 mm yields 1x1 cone input to bp
        % ecc = 30 mm yields 3x3 cone input to bp
        % Support formula extrapolated from data in Dacey ... Lee, 1999 @JRG to insert
        
        support = max(minSupport,floor(1 + (2/10)*(eccentricity)));
        
        % Standard deviation of the Gaussian for the center.  Anywhere near
        % the center the input is basically 1 cone.  Far in the periphery,
        % it will be seomthing else that we will have a function for, like
        % the support.s
        spread = 1;
        obj.sRFcenter   = fspecial('gaussian',[support,support], spread); % convolutional for now
        obj.sRFsurround = 1.3*fspecial('gaussian',[support,support], 10*spread); % convolutional for now
        
        % How do spatially sample the positions?
        obj.cellLocation = [];
        
end
