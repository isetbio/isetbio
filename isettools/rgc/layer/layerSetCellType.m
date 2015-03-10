function layerSetCellType(layer,type,theta, varargin)
% Set layer parameters to one of several predefined types
%
%   layersSetCellType(layer,type, theta, varargin)
%
% These should be cell types, such as on-parasol, small bistratified,
% off-midget, and so forth.  They could also have a parameter, such as
% eccentricity on the retina (in mm or um).
%
% Example:
%   We access this only through the rgcLayer call (addLayer)
%   rather than directly. For example, create a layer of on parasol cells:
%       rgcP = rgcParameters;
%       rgcP.addLayer('on parasol');  
%
% See also: rgcLayer, rgcParameters
%
% (c) 2010 Stanford Synapse Team
%
% TODO: 
%   Consider specifying parameter/value pairs in varargin, such as
%   eccentricity. 
%   More realistic parameter values. (e.g., cell spacing, RF size, cone
%   weights.) Should these depend on cone spacing? Eccentricty? 


if (notDefined('layer') || ~isequal(class(layer),'rgcLayer'))
    error('object of the rgcLayer class required.');
end
if notDefined('type'),  type = 'default'; end
if notDefined('theta'), theta = 0;        end

type = stringFormat(type);

coneSpacing = layer.parent.get('coneSpacing');

% get the relative frequencies of each cone type. we will need this to
% scale the coneWeights. Perhaps there should be a get function to do this.
ct = sensorGet(layer.parent.sensor, 'cone type');
ct = ct(:);
lmsFraction = [sum(ct == 2) sum(ct == 3) sum(ct == 4)]/numel(ct);	% should be determined by a get function

% Probably broke many things
% The parameters are not set well inside this switch statement.  Nearly
% everything is wrong.  (BW).
% TODO:
%  We need to make the rfCoefs specific to the cone types.
%  Find an easy way to adjust the sizes?  Maybe that's easy already.
switch(type)
    case {'default','onmidget'}
        layer.set('name',strcat(layer.get('name'),'on midget'));
        
        % One to one between coneSpacing and cellSpacing:
        layer.set('cellSpacing',layer.parent.get('coneSpacing'));
        rgcFindGoodDistanceParameters(layer,'yes');
        
        % Each cell type  should set rfCoefs and sd
        rfCoefs        = [1 -1];   % weights
        sdC            = 2;        % um 
        sdC =  (sdC/coneSpacing);  % Change to cone counts
        sdS =   2*sdC;             % sd of RF surround, in um?
        
        % Center and surround standard deviation and rotation
        sd{1} = [sdC, sdC, theta];  % x,y, rotation
        sd{2} = [sdS, sdS, theta];
        
        % No S-cones (EJ with G Field)
        centerWeight   = [1 1 0]/2; % this should sum to 1
        surroundWeight = [1 1 0]/2; 
        coneWeights = [centerWeight./ lmsFraction; ...
            surroundWeight./ lmsFraction];
        
    case {'offmidget'}
        layer.set('name',strcat(layer.get('name'),'off midget'));

        % One to one between coneSpacing and cellSpacing:
        layer.set('cellSpacing',layer.parent.get('coneSpacing'));

        % Connection matrix weightings.
        rgcFindGoodDistanceParameters(layer,'yes');
        
        % Each cell type  should set rfCoefs and sd
        rfCoefs        = [-1 1];   % weights
        sdC            = 2;        % um 
        sdC =  (sdC/coneSpacing);  % Change to cone counts
        sdS =   2*sdC;             % sd of RF surround, in um?
        
        % Center and surround standard deviation and rotation
        sd{1} = [sdC, sdC, theta];  % x,y, rotation
        sd{2} = [sdS, sdS, theta];
        
        % According to EJ, this class gets the S-cone
        centerWeight   = [1 1 1]/3;
        surroundWeight = [1 1 1]/3;
        coneWeights = [centerWeight./ lmsFraction; ...
            surroundWeight./ lmsFraction];

        
    case {'onparasol'}
        layer.set('name',strcat(layer.get('name'),'on parasol'));
        
        % cellSpacing 3 times bigger that coneSpacing
        layer.set('cellSpacing',coneSpacing*3);
        
        rgcFindGoodDistanceParameters(layer,'yes');

        % Each cell type  should set rfCoefs and sd
        rfCoefs        = [1 -1];   % weights
        sdC            = 6;        % um 
        sdC =  (sdC/coneSpacing);  % Change to cone counts
        sdS =   2*sdC;             % sd of RF surround, in um?
        
        % Center and surround standard deviation and rotation
        sd{1} = [sdC, sdC, theta];  % x,y, rotation
        sd{2} = [sdS, sdS, theta];
        
        % No S-cones (EJ with G Field)
        centerWeight   = [1 1 0]/2; % this should sum to 1
        surroundWeight = [1 1 0]/2; 
        coneWeights = [centerWeight./ lmsFraction; ...
            surroundWeight./ lmsFraction];


    case {'offparasol'}
        % We set the rf size to be 20% smaller in off center RGCs compared
        % to on-center (Functional Asymmetries in ON and OFF Ganglion Cells
        % of Primate Retina, Chichilnisky and Kalmar, JNS 2002;
        % http://neuro.cjb.net/content/22/7/2737.full)

        layer.set('name',strcat(layer.get('name'),'off parasol'));
        
        % cellSpacing 3 times bigger that coneSpacing
        layer.set('cellSpacing',coneSpacing*3);
        rgcFindGoodDistanceParameters(layer,'yes');
        
        % Each cell type  should set rfCoefs and sd.  The rfCoefs should
        % become a 3x2 for [L,M,S center; L,M,S surround];
        rfCoefs        = [-1 1];   % Weights   
        sdC            = 6;        % um        
        sdC =  (sdC/coneSpacing);  % Change to cone counts
        sdS =   2*sdC;             % sd of RF surround, in um?
        
        % Center and surround standard deviation and rotation. 
        sd{1} = [sdC, sdC, theta];  % x,y, rotation
        sd{2} = [sdS, sdS, theta];
        
        % No S-cones (EJ with G Field)
        centerWeight   = [1 1 0]/2; % this should sum to 1
        surroundWeight = [1 1 0]/2; 
        coneWeights = [centerWeight./ lmsFraction; ...
            surroundWeight./ lmsFraction];
        
    case {'smallbistratified'}
        layer.set('name',strcat(layer.get('name'),'small bistratified'));
        
        % Is 'surround' for small bistratified co-extensive with center, or
        % larger? see... 
        %
        % Parallel ON and OFF Cone Bipolar Inputs Establish Spatially
        % Coextensive Receptive Field Structure of Blue-Yellow Ganglion
        % Cells in Primate Retina        
        % Joanna D. Crook, Christopher M. Davenport, Beth B. Peterson, Orin
        % S. Packer, Peter B. Detwiler, and Dennis M. Dacey
        % http://www.jneurosci.org/content/29/26/8372.short
        
        % Low resolution : cellSpacing 4 times bigger that coneSpacing
        layer.set('cellSpacing',coneSpacing*4);
        rgcFindGoodDistanceParameters(layer,'yes');
        
        % Doesn't account for S-cones or anything.  Make it better.
        rfCoefs        = [-1 1];   % Weights
        sdC            = 12;       % um
        sdC =  (sdC/coneSpacing);  % Change to cone counts
        sdS =   2*sdC;             % sd of RF surround, in um?
        
        % Center and surround standard deviation and rotation.
        sd{1} = [sdC, sdC, theta];  % x,y, rotation
        sd{2} = [sdS, sdS, theta];
               
        centerWeight   = [0 0 1];
        surroundWeight = [1 1 0]/2;
        coneWeights = [centerWeight./ lmsFraction; ...
            surroundWeight./ lmsFraction];
        
        %     case {5}  % lr
        %         layer.set('name',strcat(layer.get('name'),'lr'));
        %         % Low resolution : coneSpacing 3 times bigger
        %         coneSpacing = layer.parent.get('coneSpacing');
        %         hrSpacing = coneSpacing*3;
        %         layer.set('cellSpacing',hrSpacing);
        %         rgcFindGoodDistanceParameters(layer,'yes');
        %
        %     case {6}  %hr
        %         layer.set('name',strcat(layer.get('name'),'hr'));
        %         % High resolution : coneSpacing 3 times smaller
        %         coneSpacing = layer.parent.get('coneSpacing');
        %         hrSpacing = coneSpacing/3;
        %         layer.set('cellSpacing',hrSpacing);
        %         rgcFindGoodDistanceParameters(layer,'yes');
        
    otherwise
        error('Unknown type %s\n',type);
        
end


layer.set('rf coefs',rfCoefs);
layer.set('cone weights',coneWeights);
layer.set('RF std',sd);

layer.hasCoupling = 1;
layer.hasFeedback = 1;

return;

%% Old code

% layer.set('name',strcat(layer.get('name'),'rotation'));
% % rotation
% theta = theta;
% R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
% RFcov = layer.get('RFcov');
% nC = length(RFcov);
% for jj = 1:nC
%     RFcov{jj} = (R')*RFcov{jj}*R;
% end
% layer.set('RFcov',RFcov);
% 
% 
% %
% 
% layer.set('name',strcat(layer.get('name'),'cornerdetector'));
% % corner detector
% rfCov = cell(1,2);
% rfCov{1} = rgcCovarianceMatrix(sqrt(5),  sqrt(5),  0);
% rfCov{2} = rgcCovarianceMatrix(sqrt(30), sqrt(30), 0);
% layer.set('rf cov',rfCov);
% 
% layer.set('RFsCoeffs',[1 -0.99]);
% layer.set('ActThreshParam',500);
