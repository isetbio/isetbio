function roiData = vcGetROIData(obj, roiLocs, dataType)
% Retrieve data within the region of interest (ROI)
%
% Syntax:
%   roiData = vcGetROIData(obj, roiLocs, dataType)
%
% Description:
%    Region of interest data (roiData) are retrieved from the window of any
%    of the ISET objects, (SCENE, OPTICALIMAGE, ISA or VCIMAGE). The
%    locations, roiLocs, can be selected using the routine vcROISelect.
%    These are an Nx2 matrix of rows and columns.
%
%    When roiLocs is sent in as a 4D row vector, the meaning is [col, row,
%    width, height];  The first data point is drawn from the position [col,
%    row]. For unfortunate historical reasons, which could be fought here,
%    the spatial size of the returned data are width + 1 and height + 1.
%    Thus, [col, row, 1, 1] returns four positions.
%
%    A variety of data types and windows can be chosen. These are listed
%    below under the dataTypes variable in the inputs section.
%
%    The data are returned in a matrix (XW format), roiData. The rows
%    represent each of the positions and the data across the columns
%    represent the values. The number of columns depends on whether the
%    data are wavelength, rgb, or electrons.
%
%    The code below contains examples of function usage. To access, type
%    'edit vcGetROIData.m' into the Command Window.
%
% Inputs:
%    obj      - Object. One of the supported objects. The objects include:
%               scene, opticalImage, and vcimage.
%    roiLocs  - Matrix. Matrix containing the region of interest locations.
%    dataType - String. String indicating the type of data to return, this
%               is dependant on the object. Supported strings, and their
%               defaults are broken up by their object types below:
%        scene:         photons (default) or energy
%                       'illuminant photons' or 'illuminant energy'
%                       reflectance
%        opticalimage:  photons (default) or energy
%        vcimage:       results (default) or input
%
% Outputs:
%    roiData  - Matrix. Matrix containig the desired RoI data.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    vcROISelect, sensorGET(isa, 'roielectrons');
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    05/11/18  jnm  Formatting

% Examples:
%{
    % ETTBSkip - Skipping broken examples, that when fixed, would require
    % user input...
    vci = vcGetObject('vci');
    roiLocs = vcROISelect(vci);
    result = vcGetROIData(vci, roiLocs, 'result');                 % Nx3

    scene = vcGetObject('scene');
    roiLocs = vcROISelect(scene);
    photons = vcGetROIData(scene, roiLocs, 'photons');             % Nx31
    photons = vcGetROIData(scene, roiLocs, 'illuminant photons');  % Nx31

    % Does this work for spatial-spectral, too?
    reflectance = vcGetROIData(scene, roiLocs, 'reflectance');
%}

if notDefined('obj'),  error('Object must be defined'); end

% The variable might be locs or a rect. We check and make sure it is locs.
if notDefined('roiLocs'), error('ROI locations must be defined');
elseif size(roiLocs, 2) == 4, roiLocs = ieRoi2Locs(roiLocs);
end

objType = vcGetObjectType(obj);
switch lower(objType)
    case {'scene'}
        % Read the ROI for the radiance data or the illuminant data.
        if notDefined('dataType'), dataType = 'photons'; end

        % Handle getting ROI for illuminant as well as radiance data.
        sz = sceneGet(obj, 'size');
        r = sz(1);
        c = sz(2);
        dataType = ieParamFormat(dataType);

        % Get data into XW format
        switch dataType
            case {'illuminantphotons'}
                img = sceneGet(obj, 'illuminant photons');

                % Handle spatial-spectral or not
                if isvector(img), img = repmat(img(:)', prod(sz), 1);
                else, [img, r, c] = RGB2XWFormat(img);
                end
            case {'illuminantenergy'}
                img = sceneGet(obj, 'illuminant energy');

                % Handle spatial-spectral or not
                if isvector(img), img = repmat(img(:)', prod(sz), 1);
                else, [img, r, c] = RGB2XWFormat(img);
                end
            case {'photons', 'energy'}
                img = sceneGet(obj, dataType);
                [img, r, c] = RGB2XWFormat(img);
                if isempty(img)
                    errordlg('No radiance data in scene.');
                end
            case {'luminance', 'luminanceroi'}
                img = sceneGet(obj, 'luminance');
                [img, r, c] = RGB2XWFormat(img);
                if isempty(img)
                    errordlg('No luminance data in scene.');
                end
        end

        % Should we keep the data in bounds?
        roiLocs(:, 1) = ieClip(roiLocs(:, 1), 1, r);
        roiLocs(:, 2) = ieClip(roiLocs(:, 2), 1, c);

        imgLocs = sub2ind([r, c], roiLocs(:, 1), roiLocs(:, 2));
        roiData = img(imgLocs, :);

    case {'opticalimage', 'oi'}
        if notDefined('dataType'), dataType = 'photons'; end

        data = oiGet(obj, dataType);
        if isempty(data)
            if strcmp(dataType , 'energy')
                photons = oiGet(obj, 'photons');
                wavelength = oiGet(obj, 'wavelength');
                data = Quanta2Energy(wavelength, photons);
            else
                errordlg('No data to plot.');
            end
        end

        [img, r, c] = RGB2XWFormat(data);

        % Should we keep the data in bounds?
        roiLocs(:, 1) = ieClip(roiLocs(:, 1), 1, r);
        roiLocs(:, 2) = ieClip(roiLocs(:, 2), 1, c);

        imgLocs = sub2ind([r, c], roiLocs(:, 1), roiLocs(:, 2));
        roiData = img(imgLocs, :);

    case 'vcimage'
        if notDefined('dataType'), dataType = 'results'; end

        data = imageGet(obj, dataType);  %input, result
        [img, r, c] = RGB2XWFormat(data);

        % Should we keep the data in bounds?
        roiLocs(:, 1) = ieClip(roiLocs(:, 1), 1, r);
        roiLocs(:, 2) = ieClip(roiLocs(:, 2), 1, c);

        imgLocs = sub2ind([r, c], roiLocs(:, 1), roiLocs(:, 2));
        roiData = img(imgLocs, :);
end

end