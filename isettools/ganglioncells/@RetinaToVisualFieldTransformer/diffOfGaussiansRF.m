function RF2D = diffOfGaussiansRF(p, spatialSupportDegs) 
    % Spatial support mesh
    [Xdegs,Ydegs] = meshgrid(spatialSupportDegs(:,1), spatialSupportDegs(:,2));
    Rdegs = sqrt(Xdegs.^2+Ydegs.^2);

    % Retrieve center params
    RcDegs = p(2); 
    p = real(p);
    Kc = p(1); 
    if (numel(p) > 4)
        flatTopGaussianExponent = p(5);
        if (numel(p) == 6)
            rotationDegs = p(6);
        else
            rotationDegs = 0;
        end
    else
        flatTopGaussianExponent = 1;
        rotationDegs = 0;
    end


    if (~isreal(RcDegs))
        % Elongated RF center
        RcDegsMajor = max([real(RcDegs) imag(RcDegs)]);
        RcDegs = min([real(RcDegs) imag(RcDegs)]);
    else
        % Circular RF center
        RcDegsMajor = RcDegs;
    end

    % Compute surround params
    RcDegsEllispoid = sqrt(RcDegs^2+RcDegsMajor^2);
    RsDegs = RcDegsEllispoid * p(3);
    Ks = Kc * p(4)/((RsDegs/RcDegsEllispoid)^2);

    rcToPixelResolutionRatio = RcDegs/(spatialSupportDegs(2,1)-spatialSupportDegs(1,1));
    if (rcToPixelResolutionRatio < 0.8)
       error('RcDegs is too small compared to PSF resolution');
    end 

    % Compute center RF
    %centerRF = Kc * exp(-(Rdegs/RcDegs).^2);
    centerRF = Kc * (exp(-(Xdegs/RcDegs).^2) .* exp(-(Ydegs/RcDegsMajor).^2)).^flatTopGaussianExponent;
    centerRF = imrotate(centerRF, rotationDegs+90, "bilinear", "crop");

    % Compute surround RF
    surroundRF = Ks * exp(-(Rdegs/RsDegs).^2);

    % Composite RF
    RF2D =  centerRF - surroundRF;
end