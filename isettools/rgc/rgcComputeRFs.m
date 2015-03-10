function RF = rgcComputeRFs(layer)
% Calculates the receptive field components RF(:,:,nComponent)
%
%  RF = rgcComputeRFs(layer)
%
% rgcComputeLayers reads the layer parameters and returns a matrix of the
% RF for each of the RF components in this layer.  
%
% For more information on any layer types, consult rgcAddLayerParameters,
% which creates parameters for a default layer
%
% Examples:
%   rgcP = rgcParameters;   %initializing param
%   rgcP.addLayer();
%   RF = rgcComputeRFs(rgcP.get('layer',1));
%
%   vcNewGraphWin; mesh(RF(:,:,1));
%
% (c) 2010 Stanford Synapse Team 

%% Check input
if notDefined('layer'), error('Layer structure required.'); end;

%% Create spatial RF
cov   = layer.get('RFcov');               % covariance matrix of each component (um^2)
wgts  = layer.get('RFsCoeffs');           % Center/surround weights
rfSupport = layer.get('rf grid');      % Spatial support (um)

nComponents = length(wgts);   % I think this is really the number of components

for ii = 1 : nComponents   % Center or surround, I think
    % Create gaussian and store spacing here set to coneSpacing, it has to
    % be that for convolution
    newGaussian = getGaussian(rfSupport, cov{ii});
    if ii==1
        % allocating space
        s1 = size(newGaussian,1);
        s2 = size(newGaussian,2);
        RF = zeros(s1,s2,nComponents);
        if (s1<2 || s2<2)
            error('The s parameter is too small for the %d-th layer.',ii);
        end
    end
    % Scale and store gaussian
    RF(:, :, ii) = wgts(ii) * newGaussian;
    % vcNewGraphWin; ii = 2; mesh(RF(:,:,ii))
end
% vcNewGraphWin; mesh(RF(:,:,1)); mesh(RF(:,:,2)); mesh(sum(RF,3))

% cutting it to a subpart of it, so we don't get too big a RF
% assuming 2 components
% sRF = abs(RF(:,:,1)+RF(:,:,2));

% threshold = max(sRF(:))*1e-4;
% mi = minimumIndex(sRF,threshold);
% mj = minimumIndex(sRF',threshold);
% 
% % assuming somehow symetric
% RF1 = RF(mi:end-mi+1,mj:end-mj+1,1);
% RF2 = RF(mi:end-mi+1,mj:end-mj+1,2);
% 
% RFs{1} = cat(3,RF1,RF2);
end

% What is this?  Can't it be better?
% function mi = minimumIndex(A,threshold)
%     A = A > threshold;
%     [a b] = max(A);
%     mi = min(b(b>1));
% end