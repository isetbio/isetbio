function RGCRFposMicrons = initializeWithPerfectHexLattice(theInputConeMosaic, coneToRGCDensityRatio)
% Initialize RGC RF positions with a perfect hex lattice and a desired coneToRGCDensityRatio
%
% Syntax:
%   RGCRFposMicrons = RGCRFconnector.initializeWithPerfectHexLattice(theInputConeMosaic, ...
%                                                     coneToRGCDensityRatio)
%
% Description:
%   Initialize RGC RF positions with a perfect hex lattice 
%   and a desired coneToRGCDensityRatio
%
% Inputs:
%    theInputConeMosaic     - The cone mosaic that will provide input
%    coneToRGCDensityRatio  - The desired cone to RGC density ratio
%
% Outputs:
%    RGCRFposMicrons        - The positions (in microns) of the RGC RFs 
%
% Optional key/value pairs
%   none
%   
% History:
%   5/11/2022       NPC     Wrote it
%

    % Generate RGCRF positions that lie within the limits of the input cone mosaic
    allConePositions = theInputConeMosaic.coneRFpositionsMicrons;
    minConePosX = min(allConePositions(:,1));
    minConePosY = min(allConePositions(:,2));
    maxConePosX = max(allConePositions(:,1));
    maxConePosY = max(allConePositions(:,2));

    centerMicrons = [mean(allConePositions(:,1)) mean(allConePositions(:,2))];
    sizeMicrons = max([maxConePosX-minConePosX maxConePosY-minConePosY]);

    % Generate perfect hex lattice of RGC RF positions corresponding to the
    % coneToRGCdensity ratio
    mRGCseparation = min(theInputConeMosaic.coneRFspacingsMicrons) * sqrt(2.0./(sqrt(3.0)/coneToRGCDensityRatio));
    lambdaMicrons = 0.97*mRGCseparation;
    RGCRFposMicrons = generateHexLattice(centerMicrons, sizeMicrons/2, lambdaMicrons, 'rectangular');

    idx = find(...
        (RGCRFposMicrons(:,1) >= minConePosX) & ...
        (RGCRFposMicrons(:,1) <= maxConePosX) & ...
        (RGCRFposMicrons(:,2) >= minConePosY) & ...
        (RGCRFposMicrons(:,2) <= maxConePosY));
    RGCRFposMicrons = RGCRFposMicrons(idx,:);

end

function XY = generateHexLattice(center,radius, lambda, domain) 
    upSampleFactor = 1000;
    radius = radius * upSampleFactor;
    lambda = lambda * upSampleFactor;
    circlesNum = 8+max([1 ceil(round(radius/lambda))]);

    X = []; Y = [];
    for iCircle = 1:circlesNum
        x = zeros(iCircle*6,1);
        y = zeros(iCircle*6,1);

        x(1:6) = lambda * iCircle * cos(2*pi/6.*(0:5));
        y(1:6) = lambda * iCircle * sin(2*pi/6.*(0:5));
        if iCircle > 1
            for q = 1:iCircle-1
               x0 = x(2) - q*lambda;
               radi0 = sqrt(x0^2+y(2)^2);
               theta0 = 1/3*pi + pi * 1/3 * 1/(iCircle) * q;
               x(q*6 + 1:(q*6+6)) = radi0*cos(theta0 + pi/3 .* (1:6));
               y(q*6 + 1:(q*6+6)) = radi0*sin(theta0 + pi/3 .* (1:6));
            end
        end
        X = [X; x];
        Y = [Y; y];
    end

    % Add RFposition at (0,0)
    X = [0; X]; Y = [0; Y];

    switch (domain)
        case 'circular'
            d = sqrt(sum(X.^2 + Y.^2, 2));
            idx = find(d <= radius);
        case 'rectangular'
            idx = find(...
                (X >= -radius) & ...
                (X <=  radius) & ...
                (Y >= -radius) & ...
                (Y <=  radius));
    end

    X = X(idx)/upSampleFactor;
    Y = Y(idx)/upSampleFactor;
    XY = [center(1)+X(:) center(2)+Y(:)];
end


