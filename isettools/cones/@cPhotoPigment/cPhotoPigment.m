classdef cPhotoPigment < receptorPigment
% Class for single cone photopigment and related properties
%
% Syntax:
%	pigment = cPhotoPigment;
%
% Description:
%    This class contains properties for the photopigment absorption
%    properties of a single cone cell. This class is derived from the 
%    @photoPigment, which handles all the spectral properties. The main
%    function of the @cPhotoPigment class is to handle geometry for a
%    disk-shaped cone aperture. It is to be used with the (new) @cMosaic class.
%
%    For the full cone mosaic, see the @cMosaic class
%
%    See extensive header comments in the @receptorPigment class for more
%    information.
%
% Input:
%	 None required.
%
% Output:
%    pigment          - The created cPhotoPigment object.
%   
% Optional key/value pairs:
%	 'wave'           - Vector of wavelengths in nm (400:10:31).
%    'opticalDensity' - Three vector of optical densities for L, M and
%                       S cone photopigment. Default: [0.5 0.5 0.4].
%    'absorbance'     - L, M and S cone absorbance spectra. Default
%                       empty, in which case these are obtained through
%                       routine coneAbsorbanceReadData.
%    'peakEfficiency' - Quantal efficiency for isomerizations for
%                       L, M and S cones. Default [2 2 2]/3.
%    'diameter'       - Diameter of the disk-shaped light collecting aperture. Default 4.0/sqrt(pi))*1e-6.
%
% See Also:
%   photoPigment, cMosaic, Macular, Lens
%

% History:
%    12/07/20  NPC Wrote it

properties  % public properties related to the geometry of the cone aperture
    % diameter of the light collecting cone aperture in meters
    diameter;
end


properties (Dependent)
    % area - The area of the light collecting aperture.
    area;
end


methods  % public methods
    % constructor
    function obj = cPhotoPigment(varargin)
        p = inputParser;
        p.KeepUnmatched = true;
        p.addParameter('wave', 400:10:700, @isnumeric);
        p.addParameter('opticalDensity', [0.5 0.5 0.4], @isnumeric);
        p.addParameter('absorbance', [], @isnumeric);
        p.addParameter('peakEfficiency', [2 2 2]/3, @isnumeric);
        p.addParameter('diameter', (4.0/sqrt(pi))*1e-6, @isnumeric);
        p.parse(varargin{:});

        % Call the super-class constructor.
        obj = obj@receptorPigment(...
            'wave', p.Results.wave(:), ...
            'opticalDensity', p.Results.opticalDensity, ...
            'absorbance', p.Results.absorbance, ...
            'peakEfficiency', p.Results.peakEfficiency);
        
        % set object properties
        obj.diameter = p.Results.diameter;
        
        % Assert that property dimensions are consistent
        assert(numel(obj.opticalDensity) == numel(obj.peakEfficiency), ...
            sprintf('optical density dimensionality does not match that of peak efficiency'));
        
        % Assert that property dimensions are consistent
        if (~isempty(p.Results.absorbance))
            assert(numel(obj.opticalDensity) == size(obj.absorbance,2), ...
                sprintf('optical density dimensionality does not match that of absorbance'));  
        end
        
    end


    function val = get.area(obj)
        % Compute  the light-collecting area 
        %
        % Syntax:
        %   obj = get.area(obj)
        %
        % Description:
        %    Compute  the light-collecting area 
        %
        % Inputs:
        %    obj - The cPhotoPigment object
        %
        % Outputs:
        %    val - The area value for obj
        %
        % Optional key/value pairs:
        %    None.
        %
        val = pi * (0.5*obj.diameter)^2;
    end

end

methods (Static)
    % When we have them, they go here
end
end