function RF2D = diffOfGaussiansRF(p, spatialSupportDegs) 
    % Spatial support mesh
    [Xdegs,Ydegs] = meshgrid(spatialSupportDegs(:,1), spatialSupportDegs(:,2));
    Rdegs = sqrt(Xdegs.^2+Ydegs.^2);

    % Retrieve center params
    Kc = p(1); 
    RcDegs = p(2); 
    
    % Retrieve surround params
    RsDegs = RcDegs * p(3);
    Ks = Kc * p(4)/((RsDegs/RcDegs)^2);

    rcToPixelResolutionRatio = RcDegs/(spatialSupportDegs(2,1)-spatialSupportDegs(1,1));
    if (rcToPixelResolutionRatio < 0.8)
       error('RcDegs is too small compared to PSF resolution');
    end 

    % Compute center RF
    centerRF = Kc * exp(-(Rdegs/RcDegs).^2);

    % Compute surround RF
    surroundRF = Ks * exp(-(Rdegs/RsDegs).^2);

    % Composite RF
    RF2D =  centerRF - surroundRF;
end