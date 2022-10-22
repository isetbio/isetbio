function eccDegs = temporalEquivalentEccentricityFromEccentricity(obj, varargin)
    p = inputParser;
    p.addParameter('forEccDegsMatrix', [], @(x)(isempty(x) || isnumeric(x)));
    p.parse(varargin{:});
    eccDegs = p.Results.forEccDegsMatrix;

    if (isempty(eccDegs))
        eccDegs = obj.eccentricityDegs;
    elseif (size(eccDegs,2) ~= 2)
        error('eccDegs must be an N x 2 matrix of (x,y) eccentricities.');
    end
    dataPoints = size(eccDegs,1);

    for iPoint = 1:dataPoints
        % The vertical ecc
        verticalEcc = eccDegs(iPoint,2);

        % The horizontal ecc
        if strcmp(obj.horizontalRetinalMeridian, RGCmodels.Watson.constants.nasalMeridian)
            % On the nasal meridian
            switch (obj.temporalEquivantEccentricityFactor)
                case 'WatanabeRodieckBased'
                    % Nasal meridian: multiply by 0.61 - See Watanabe & Rodieck 1989
                    horizontalEcc = 0.61 * eccDegs(iPoint,1);
                case 'ISETBioMosaicsBased'
                    % Nasal meridian: multiply by 0.75 - to match our asymmetry
                    horizontalEcc = 0.72 * eccDegs(iPoint,1);
            end
        else
            % We are already on the temporal meridian
            horizontalEcc = eccDegs(iPoint,1);
        end

        % The radial temporal equivalent eccentricity
        radialEquivalentEccentricity = sqrt(horizontalEcc^2 + verticalEcc^2);

        % Compute the vector temporal equivalent eccentricity
        theta = angle(eccDegs(iPoint,1) + 1j*eccDegs(iPoint,2));
        eccDegs(iPoint,:) = radialEquivalentEccentricity  * [-cos(theta) sin(theta)];
    end % iPoint
end