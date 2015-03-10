%% s_MicroLensPhaseSpace
%
% Under development
%
%    Illustration of the steps used in the phase-space calculatios for
%    Micro-lens analysis.  The computations use Wigner phase space (without
%    diffraction).  The individual routines and optical properties are set
%    by hand in this script.
%
% Peter Catrysse to comment and update.
%

% Initialization Parameters
nAir = 1;                       % Index of refraction in air
widthPS = 6;                    % [um], Wigner Phase-Space size in microns

% Micro-Lens Parameters
f = 9;                          % [um]
lensOffset = 3.23; % 0, 0.82, 1.65, 2.45, 3.23

% Pixel Parameters
pixelWidth = 3.3;               % microns

% Source Parameters
lambda = 0.5;		            % [um]
sourceWidth = pixelWidth;       % [um]   Why is the source width = pixel width?
sourcefNumber = 2.8;
sourceChiefRayAngle = 20; % 0, 5, 10, 15, 20

% Wigner PS grids
[X,P] = mlCoordinates(-widthPS,widthPS,nAir,lambda,'angle');
x = X(1,:); p = P(:,1);

% Lambertian source
sourceNA = sin(atan(1/(2*sourcefNumber)));
[W,X,P] = mlSource(-sourceWidth/2,sourceWidth/2, ...
    sin(sourceChiefRayAngle/180*pi) - sourceNA, ...
    sin(sourceChiefRayAngle/180*pi) + sourceNA,X,P,'angle');

figure(1); subplot(3,3,2); imagesc(X(1,:), P(:,1),W)
set(gca,'YDir','normal');
colormap('hot'); title('Source'); xlabel('x (\mum)'); ylabel('n sin(\theta)');
figure(1); subplot(3,3,3); plot(sum(W,2)/max(sum(W,2)),p); grid;
title('Angular Intensity'); ylabel('n sin(\theta)');
figure(1); subplot(3,3,1); plot(x,sum(W,1)/max(sum(W,1))); grid;
title('Positional Intensity'); xlabel('x (\mum)');

% Lens
[V_,X,P] = mlLens(f,lambda,W,X,P,'non-paraxial','angle');

figure(1); subplot(3,3,5); imagesc(X(1,:),P(:,1),V_); set(gca,'YDir','normal');
% mesh(X,P,V_) 
colormap('hot'); title('After lens'); xlabel('x (\mum)'); ylabel('n sin(\theta)');
figure(1); subplot(3,3,6); plot(sum(V_,2)/max(sum(V_,2)),p); grid;
title('Angular Intensity'); ylabel('n sin(\theta)');
figure(1); subplot(3,3,4); plot(x,sum(V_,1)/max(sum(V_,1))); grid;
title('Positional Intensity'); xlabel('x (\mum)');

% Lens displacement
[V__,X,P] = mlDisplacement(lensOffset,V_,X,P,'non-paraxial','angle');

figure(1); subplot(3,3,8); imagesc(X(1,:),P(:,1),V__); set(gca,'YDir','normal');
% surf(X,P,V_) 
colormap('hot'); title('After lens'); xlabel('x (\mum)'); ylabel('n sin(\theta)');
figure(1); subplot(3,3,9); plot(sum(V__,2)/max(sum(V__,2)),p); grid;
title('Angular Intensity'); ylabel('n sin(\theta)');
figure(1); subplot(3,3,7); plot(x,sum(V__,1)/max(sum(V__,1))); grid;
title('Positional Intensity'); xlabel('x (\mum)');

% Propagation over distance f
[W_,X,P] = mlPropagate(f,lambda,V__,X,P,'non-paraxial','angle');

figure(1); subplot(3,3,8); imagesc(X(1,:),P(:,1),W_); set(gca,'YDir','normal');
colormap('hot'); title('After distance f'); xlabel('x (\mum)'); ylabel('n sin(\theta)');
figure(1); subplot(3,3,9); plot(sum(W_,2)/max(sum(W_,2)),p); grid;
title('Angular Intensity'); ylabel('n sin(\theta)');
figure(1); subplot(3,3,7); plot(x,sum(W_,1)/max(sum(W_,1))); grid;
title('Positional Intensity'); xlabel('x (\mum)');

% Optical Efficiency
IrradianceIn = sum(W,1); etendueIn = sum(IrradianceIn(find(abs(x)<(sourceWidth/2))));
IrradianceOut = sum(W_,1); etendueOut = sum(IrradianceOut(find(abs(x)<(sourceWidth/2))));
E = etendueOut/etendueIn;

% Map of Irradiance
pixelIrradiance = repmat(sum(W_,1)/max(sum(W_,1)),length(x),1).*rot90(repmat(sum(W_,1)/max(sum(W_,1)),length(x),1));
pixelIrradiance(find((abs(x) < 1.7)&(abs(x) > 1.6)),find((abs(x) < 1.7)&(abs(x) > 1.6))) = 1;
pixelIrradiance(find((abs(x) < 1.7)&(abs(x) > 1.6)),:) = 1;
pixelIrradiance(:,find((abs(x) < 1.7)&(abs(x) > 1.6))) = 1;
figure(2); imagesc(x,x,pixelIrradiance); colormap('hot');
