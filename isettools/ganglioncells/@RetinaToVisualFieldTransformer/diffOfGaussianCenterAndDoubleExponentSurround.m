function RF2D = diffOfGaussianCenterAndDoubleExponentSurround(p, spatialSupportDegs)

    % Spatial support mesh
    [Xdegs,Ydegs] = meshgrid(spatialSupportDegs(:,1), spatialSupportDegs(:,2));
    Rdegs = sqrt(Xdegs.^2+Ydegs.^2);

    % Retrieve center params
    Kc = p(1); 
    RcDegs = p(2); 
    
    % Retrieve surround params: sum of 2 exponentials (wide field + narrow field)
    KsToKcPeakRatio = p(3);         % range [0.01 to 1]
    narrowToWideVolumeRatio = p(4); % range [0.2 to 1.0]
    RwideDegs = p(5);               % radius at which wide sensitivity drops to 10%
    RnarrowToRwideRatio = p(6);     % range [0.01 to 1.0]
    
    RnarrowDegs = RwideDegs * RnarrowToRwideRatio;

    %Knarrow = 1.0;
    %narrowToWideVolumeRatio = Knarrow/Kwide * (RnarrowDegs/RwideDegs);
    Kwide = (RnarrowDegs/RwideDegs) / narrowToWideVolumeRatio;

    % KsToKcPeakRatio = Ks*(1+Kwide)/Kc;
    Ks = (KsToKcPeakRatio * Kc)/(1+Kwide);

    rcToPixelResolutionRatio = RcDegs/(spatialSupportDegs(2,1)-spatialSupportDegs(1,1));
    if (rcToPixelResolutionRatio < 0.8)
       error('RcDegs is too small compared to PSF resolution');
    end 

    % Compute center RF
    centerRF = Kc * exp(-(Rdegs/RcDegs).^2);

    % Compute surround RF (double exponential)
    surroundRF = Ks * ( exp(-2.3*Rdegs/RnarrowDegs) + Kwide * exp(-2.3*Rdegs/RwideDegs));

    % Composite RF
    RF2D =  centerRF - surroundRF;
end