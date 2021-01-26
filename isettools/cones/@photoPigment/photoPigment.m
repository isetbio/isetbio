classdef photoPigment < receptorPigment
% Class for single cone photopigment and related properties
%
% Syntax:
%	pigment = photoPigment;
%
% Description:
%    This class contains properties for the photopigment absorption
%    properties of a single cone cell. This class is derived from the 
%    @photoPigment, which handles all the spectral properties. The main
%    function of the @photoPigment class is to handle geometry for a
%    rectangular-shaped cone aperture. It is to be used with the (old) @coneMosaic class.
%
%    For the full cone mosaic, see the @coneMosaic and @coneMosaicHex classes.
%
%
% Input:
%	 None required.
%
% Output:
%    pigment          - The created photoPigment object.
%   
% Optional key/value pairs:
%	 'wave'           - Vector of wavelengths in nm (400:10:31).
%    'opticalDensity' - Three vector of peak optical densities for L, M and
%                       S cone photopigment. Default: [0.5 0.5 0.4].
%    'absorbance'     - L, M and S cone absorbance spectra. Default
%                       empty, in which case these are obtained through
%                       routine coneAbsorbanceReadData.
%    'peakEfficiency' - Quantal efficiency for isomerizations for
%                       L, M and S cones. Default [2 2 2]/3.
%    'width'          - Cone width (including gap between cones) in
%                       meters. Default 2e-6.
%    'height'         - Cone height (including gap between cones) in
%                       meters. Default 2e-6.
%    'pdWidth'        - Collecting area width in meters (default 2e-6)
%    'pdHeight'       - Collecting area height in meters (default 2e-6)
%
% Notes:
%    * [NOTE: DHB - Need to explain about width and height, pdWidth and
%      pdHeight and how these are used. Perhaps even simplify code not
%      to have both.]
%
% See Also:
%    t_conePhotoPigment, cPhotoPigment, coneMosaic, Macular, Lens
%

% History:
%    xx/xx/16  HJ   ISETBIO Team, 2016
%    02/15/18  jnm  Formatting
%    12/18/20  dhb  Comments.  Add quantalEfficiency property.

properties  %  % public properties related to the geometry of the cone aperture
    % width - cone width (including gap) in meters
    width;

    % height - cone height (including gap) in meters
    height;

    % pdWidth - photodetector width in meters
    pdWidth;

    % pdHeight - photodetector height in meters
    pdHeight;
end


properties (Dependent)
    % area - The area of the object. Calculated by width * height
    area;

    % pdArea - The pdArea of the object. Calculated by pdWidth * pdHeight
    pdArea;

    % gapWidth - the width of the gap. Calculated by width - pdWidth
    gapWidth;

    % gapHeight - The height of the gap. Calculated by height - pdHeight
    gapHeight;
end


methods  % public methods
    % constructor
    function obj = photoPigment(varargin)
        p = inputParser;
        p.KeepUnmatched = true;
        p.addParameter('wave', 400:10:700, @isnumeric);
        p.addParameter('opticalDensity', [0.5 0.5 0.4], @isnumeric);
        p.addParameter('absorbance', [], @isnumeric);
        p.addParameter('peakEfficiency', [2 2 2]/3, @isnumeric);
        p.addParameter('width', 2e-6, @isnumeric);
        p.addParameter('height', 2e-6, @isnumeric);
        p.addParameter('pdWidth', 2e-6, @isnumeric);
        p.addParameter('pdHeight', 2e-6, @isnumeric);
        p.parse(varargin{:});

        % Call the super-class constructor.
        obj = obj@receptorPigment(...
            'wave', p.Results.wave(:), ...
            'opticalDensity', p.Results.opticalDensity, ...
            'absorbance', p.Results.absorbance, ...
            'peakEfficiency', p.Results.peakEfficiency);
        
        % set object properties
        obj.width = p.Results.width;
        obj.height = p.Results.height;
        obj.pdWidth = p.Results.pdWidth;
        obj.pdHeight = p.Results.pdHeight;
    end

    function val = get.area(obj)
        % Retrieve photo pigment object's area
        %
        % Syntax:
        %   obj = get.area(obj)
        %
        % Description:
        %    Retrieve the area from the photoPigment object obj
        %
        % Inputs:
        %    obj - The photoPigment object
        %
        % Outputs:
        %    val - The area value for obj
        %
        % Optional key/value pairs:
        %    None.
        %
        val = obj.width * obj.height;
    end

    function val = get.gapWidth(obj)
        % Retrieve photo pigment object's gap width value
        %
        % Syntax:
        %   obj = get.gapWidth(obj)
        %
        % Description:
        %    Retrieve the gap width from the photoPigment object obj
        %
        % Inputs:
        %    obj - The photoPigment object
        %
        % Outputs:
        %    val - The gap width value for obj
        %
        % Optional key/value pairs:
        %    None.
        %
        val = obj.width - obj.pdWidth;
    end

    function val = get.gapHeight(obj)
        % Retrieve photo pigment object's gap height value
        %
        % Syntax:
        %   obj = get.gapHeight(obj)
        %
        % Description:
        %    Retrieve the gap height from the photoPigment object obj
        %
        % Inputs:
        %    obj - The photoPigment object
        %
        % Outputs:
        %    val - The gap height value for obj
        %
        % Optional key/value pairs:
        %    None.
        %
        val = obj.height - obj.pdHeight;
    end

    function val = get.pdArea(obj)
        % Retrieve photo pigment object's pd area value
        %
        % Syntax:
        %   obj = get.pdArea(obj)
        %
        % Description:
        %    Retrieve the pd area from the photoPigment object obj
        %
        % Inputs:
        %    obj - The photoPigment object
        %
        % Outputs:
        %    val - The pd area value for obj
        %
        % Optional key/value pairs:
        %    None.
        %
        val = obj.pdWidth * obj.pdHeight;
    end
end

methods (Static)
    % When we have them, they go here
end
end