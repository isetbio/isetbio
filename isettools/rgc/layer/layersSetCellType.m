function layersSetCellType(layer,type,theta)
% Set layer parameters to one of several predefined types
%
%   layersSetCellType(layer,type)
%
% These should be cell types, such as on-parasol, small bistratified,
% off-midget, and so forth.  They could also have a parameter, such as
% eccentricity on the retina (in mm or um).
%
% Example:
%   We access this only through the rgcLayer call (addLayer)
%   rather than directly.
%
% See also: rgcLayer
%
% (c) 2010 Stanford Synapse Team

if (notDefined('layer') || ~isequal(class(layer),'rgcLayer'))
    error('object of the rgcLayer class required.');
end
if notDefined('type'),  type = 'default'; end
if notDefined('theta'), theta = 0;        end

type = stringFormat(type);

coneSpacing = layer.parent.get('coneSpacing');

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
        sd{1} = [sdC, sdC, 0];  % x,y, rotation
        sd{2} = [sdS, sdS, 0];
        
    case {'offmidget'}
        layer.set('name',strcat(layer.get('name'),'off midget'));
        
        % One to one between coneSpacing and cellSpacing:
        rgcFindGoodDistanceParameters(layer,'yes');
        
        % Each cell type  should set rfCoefs and sd
        rfCoefs        = [1 -1];   % weights
        sdC            = 2;        % um 
        sdC =  (sdC/coneSpacing);  % Change to cone counts
        sdS =   2*sdC;             % sd of RF surround, in um?
        
        % Center and surround standard deviation and rotation
        sd{1} = [sdC, sdC, 0];  % x,y, rotation
        sd{2} = [sdS, sdS, 0];
        
    case {'onparasol'}
        layer.set('name',strcat(layer.get('name'),'on parasol'));
        
        % Low resolution : cellSpacing 3 times bigger that coneSpacing
        layer.set('cellSpacing',coneSpacing*3);
        rgcFindGoodDistanceParameters(layer,'yes');

        % Each cell type  should set rfCoefs and sd
        rfCoefs        = [1 -1];   % weights
        sdC            = 6;        % um 
        sdC =  (sdC/coneSpacing);  % Change to cone counts
        sdS =   2*sdC;             % sd of RF surround, in um?
        
        % Center and surround standard deviation and rotation
        sd{1} = [sdC, sdC, 0];  % x,y, rotation
        sd{2} = [sdS, sdS, 0];

    case {'offparasol'}
        % We set the rf size to be 20% smaller in off center RGCs compared
        % to on-center (Functional Asymmetries in ON and OFF Ganglion Cells
        % of Primate Retina, Chichilnisky and Kalmar, JNS 2002;
        % http://neuro.cjb.net/content/22/7/2737.full)

        layer.set('name',strcat(layer.get('name'),'off parasol'));
        
        % Low resolution : cellSpacing 3 times bigger that coneSpacing
        layer.set('cellSpacing',coneSpacing*3);
        rgcFindGoodDistanceParameters(layer,'yes');
        
        % Each cell type  should set rfCoefs and sd.  The rfCoefs should
        % become a 3x2 for [L,M,S center; L,M,S surround];
        rfCoefs        = [-1 1];   % Weights   
        sdC            = 6;        % um        
        sdC =  (sdC/coneSpacing);  % Change to cone counts
        sdS =   2*sdC;             % sd of RF surround, in um?
        
        % Center and surround standard deviation and rotation. 
        sd{1} = [sdC, sdC, 0];  % x,y, rotation
        sd{2} = [sdS, sdS, 0];
        
    case {'smallbistratified'}
        layer.set('name',strcat(layer.get('name'),'small bistratified'));
        
        % Low resolution : cellSpacing 4 times bigger that coneSpacing
        layer.set('cellSpacing',coneSpacing*4);
        rgcFindGoodDistanceParameters(layer,'yes');
        
        % Doesn't account for S-cones or anything.  Make it better.
        rfCoefs        = [-1 1];   % Weights
        sdC            = 12;        % um
        sdC =  (sdC/coneSpacing);  % Change to cone counts
        sdS =   2*sdC;             % sd of RF surround, in um?
        
        % Center and surround standard deviation and rotation.
        sd{1} = [sdC, sdC, 0];  % x,y, rotation
        sd{2} = [sdS, sdS, 0];
        
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
