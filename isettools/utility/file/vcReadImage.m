function [photons, illuminant, basis, comment, mcCOEF, linearizedImage] = vcReadImage(...
    fullname, imageType, varargin)
% Read image color data, return multispectral photons
%
% Syntax:
%   [photons, illuminant, basis, comment, mcCOEF, linearizedImage] = ...
%             vcReadImage([fullname], [imageType], [varargin])
%
% Description:
%    The image data in fullname are converted into photons. The other
%    parameters can be returned if needed. This routine is called pretty
%    much only by sceneFromFile.
%
%    There are several different image file types. This program tries to
%    determine the type from the file name. If the function is unable to
%    determine the filetype from the name, the user will be queried.
%
%    Examples are located within the code. To access the examples, type
%    'edit vcReadImage.m' into the Command Window.
%
% Inputs:
%    fullname   - (Optional) Either a file name or possible RGB data read
%                 from a file. An empty input filename produces a return, 
%                 with no error message, to work smoothly with canceling
%                 vcSelectImage. Default is to select an image based on
%                 imageType using vcSelectImage.
%    imageType  - (Optional) The type of input data. Default is 'rgb'.
%                 There are two general types: {'rgb' 'unispectral'
%                 'monochrome'} and {'multispectral', 'hyperspectral'}.
%                 These have different options for varargin.
%           {'rgb' 'unispectral' 'monochrome'}
%                      varargin{1} can be either the file name of a display
%                      (displayCreate) structure, or the display structure
%                      itself. The data in the proscribed format is
%                      returned as photons estimated by puttig the data
%                      into the display framebuffer. If there is no display
%                      calibration file, we arrange the values so that the
%                      display code returns the same RGB values as the
%                      original file.
%                      if 'rgb', then varargin{2} may contain a
%                      doSub flag. If true, the input image will be
%                      converted to a subpixel rendered image. The
%                      conversion is done by replacing each pixel with the
%                      subpixel structure of the display. Seethe function
%                      displayCompute for the details of this conversion.
%                      By default, doSub is set to False. Note that, for
%                      large image, turning on doSub could be extremely an
%                      slow process.
%           {'multispectral', 'hyperspectral'}
%                      The data are stored as coefficients and
%                      basis functions. We build the spectral
%                      representation here. These, along with a comment and
%                      measurement of the scene illuminant (usually
%                      measured using a PhotoResearch PR-650 spectral
%                      radiometer) can be returned.
%
% Outputs:
%    photons    - RGB format of photon data (r, c, w)
%    illuminant - An illuminant structure
%    basis      - A structure containing the basis functions for a
%                 multispectral SPD
%    comment    - A comment string
%    mcCOEF     - Coefficients for basis functions for multispectral SPD
%    linearizedImage - RGB image after linearization by display gamma.
%                 Only set to something interesting on some ways through
%                 this routine.  Otherwise empty.
%
% Optional key/value pairs:
%    Needs to be populated.
%
% See Also:
%    displayCompute, v_displayLUT, vcSelectImage, sceneFromFile
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    11/29/17  jnm  Formatting & notes
%    01/29/18  jnm  Formatting update to match Wiki.
%    02/14/19  dhb, lz  Return linearized image too.

% Examples:
%{
    fName = fullfile(isetbioDataPath, 'images', 'rgb', 'eagle.jpg');
    photons = vcReadImage(fName, 'rgb');
    vcNewGraphWin;
    imageSPD(photons, 400:10:700);
%}
%{
    thisDisplay = displayCreate('LCD-Apple.mat');
    fName = fullfile(isetbioDataPath, 'images', 'rgb', 'eagle.jpg');
    photons = vcReadImage(fName, 'rgb', thisDisplay);
    vcNewGraphWin;
    imageSPD(photons, displayGet(thisDisplay, 'wave'));
%}

if notDefined('imageType'), imageType = 'rgb'; end
if notDefined('fullname')
    [fullname, imageType] = vcSelectImage(imageType);
end

if isempty(fullname), photons = []; return; end

% These are loaded for a file, when they are returned.
mcCOEF = [];
comment = '';

% This is sometimes set to the linearized image
linearizedImage = [];

% Initialize other parameters to empty that are sometimes set
illuminant = [];
basis = [];

imageType = ieParamFormat(imageType);

switch lower(imageType)
    
    case {'rgb', 'unispectral', 'monochrome'}
        if isempty(varargin) || isempty(varargin{1})
            dispCal = [];
        else
            dispCal = varargin{1};
        end
        
        % doSub indicates if we should use subpixel rendering techniques
        if length(varargin) > 1
            doSub = varargin{2};
        else
            doSub = false;
        end
        
        % if we do subpixel rendering, the user could specify how many
        % samples per dixel to be used in scene generation
        if length(varargin) > 2
            sz = varargin{3};
        else
            sz = [];
        end
        
        % [Note: HJ - Programming Note: Scene rendered without subpixel
        % (doSub=false) could be different from the one rendered with
        % subpixel of size [m n], where m and n indicates the number of
        % pixels per dixel in row and columns. This kind of inconsistancy
        % will occur especially when pixels per dixel is not [1 1]. To be
        % more accurate, please turn off subpixel rendering if the samples
        % per dixel gets very small.]
        
        % Read the image data and convert them to double
        if ischar(fullname)
            inImg = double(imread(fullname));
        else
            inImg = double(fullname);
        end
        
        % A rgb image.
        if isempty(dispCal)
            if ismatrix(inImg), inImg = repmat(inImg, [1 1 3]); end
            if ndims(inImg) ~= 3
                error('Bad number of dimensions %.0f of image', ...
                    ndims(img));
            end
            % If there is no display calibration file, we arrange the
            % photon values so that the scene window shows the same RGB
            % values as in the original file.
            fprintf('[%s]: Assuming input data is 8 bit\n', mfilename);
            fprintf('[%s]: Using block matrix primaries\n', mfilename);
            [xwImg, r, c, ~] = RGB2XWFormat(inImg / 255);
            
            % Prevent DR > 10, 000. See ieCompressData.
            xwImg = ieClip(xwImg, 1e-3, 1);
            
            % When we render the RGB data in xwImg, they are multipled by
            % the colorBlockMatrix. By storing the photons this way, the
            % displayed image in the scene window will be the same as the
            % original RGB image.
            photons = xwImg * pinv(colorBlockMatrix(31));
            
        else
            % The user sent a display calibration file. If the user sent a
            % string, read the file. If the user sent in the display
            % structure, set it.
            if ischar(dispCal)
                d = displayCreate(dispCal);
            elseif isstruct(dispCal) && isequal(dispCal.type, 'display')
                d = dispCal;
            else
                error('Unknown diplay structure');
            end
            
            % Get the parameters from the display
            wave = displayGet(d, 'wave');  % Primary wavelengths
            gTable = displayGet(d, 'gamma table');
            
            np = displayGet(d, 'n primaries');
            if ismatrix(inImg), inImg = repmat(inImg, [1 1 np]); end
            
            % Pad the third dimension of inImg to nprimaries. The padding
            % is mainly for multi-primary displays, e.g. rgbw display.
            % Padding zeros might not be a good idea, but we don't have a
            % general solution
            inImg = padarray(inImg, [0 0 np-size(inImg, 3)], 0, 'post');
            assert(size(inImg, 3) == np, 'bad image size');
            
            % Check whether the gTable has enough entries for this image
            if max(inImg(:)) > size(gTable, 1)
                error('Img exceeds gTable');
            elseif max(inImg(:)) <= 1
                % DAC values are [0, 2^nBits - 1]
                inImg = round(inImg * (size(gTable, 1) - 1));
            elseif max(inImg(:)) <= 255
                % We believe this is an 8 bit image. We check whether the
                % gTable is 8 or 10 or whatever. If it is not 8 bit, then
                % we stretch the image values out to span the same range as
                % the gTable.
                s = size(gTable, 1);
                if s > 256
                    fprintf(['[%s] Assuming an 8bit image and a %d bit' ...
                        ' LUT\n'], mfilename, log2(s));
                    inImg = round(inImg / 255 * (s - 1));
                end
            end
            
            % Convert DAC values to linear intensities for the channels.
            inImg = ieLUTDigital(inImg, gTable);
            linearizedImage = inImg;
            
            % Subpixel rendering
            if doSub
                inImg = displayCompute(d, inImg, sz);
            end
            spd = displayGet(d, 'spd');   % Primary SPD in energy
            
            [xwImg, r, c] = RGB2XWFormat(inImg);
            
            % Convert energy units to quanta. This step could be slow, 
            % especially when we use sub-pixel sampling
            if numel(xwImg) < ieSessionGet('image size threshold') ...
                    || ~ieSessionGet('waitBar') % small image
                % compute directly
                photons = Energy2Quanta(wave, (xwImg * spd')')';
            else
                % now we have too many pixels, loop over wavelength and
                % print the progress
                wBar = waitbar(0, 'Computing radiance of scene');
                photons = zeros(size(xwImg, 1), length(wave));
                for ii = 1 : length(wave)
                    photons(:, ii) = Energy2Quanta(wave(ii), ...
                                     (xwImg * spd(ii, :)')')';
                    waitbar(ii / length(wave), wBar);
                end
                delete(wBar);
            end
            
            am = displayGet(d, 'ambient spd');
            amQuanta = Energy2Quanta(wave, am(:));
            photons = bsxfun(@plus, photons, amQuanta');
        end
        photons = XW2RGBFormat(photons, r, c);
        
    case {'multispectral', 'hyperspectral'}
        % These are always there. Illuminant should be there, too. But
        % sometimes it isn't, so we check below, separately.
        
        % See if the representation is a linear model with basis functions
        variables = whos('-file', fullname);
        if ieVarInFile(variables, 'mcCOEF')
            disp('Reading multispectral data with mcCOEF.')
            
            % Make this a function.
            % [photons, basis] = ieReadMultispectralCoef(fullname);
            
            % The data are stored using a linear model
            load(fullname, 'mcCOEF', 'basis', 'comment');
            
            % Resample basis functions to the user specified wavelength
            % list. vcReadImage(fullname, 'multispectral', [400:20:800]);
            if ~isempty(varargin) && ~isempty(varargin{1})
                oldWave = basis.wave; %#ok
                newWave = varargin{1};
                nBases = size(basis.basis, 2);
                extrapVal = 0;
                newBases = zeros(length(newWave), nBases);
                for ii = 1:nBases
                    newBases(:, ii) = interp1(oldWave(:), ...
                        basis.basis(:, ii), newWave(:), 'linear', ...
                        extrapVal);
                end
                basis.basis = newBases;
                basis.wave = newWave;
            end
            
            % The image data should be in units of photons
            photons = imageLinearTransform(mcCOEF, basis.basis');
            % vcNewGraphWin;
            % imageSPD(photons, basis.wave);
            
            % These lines are left in because there must be different file
            % types out there somewhere. Sometimes we stored the mean, and
            % sometimes we didn't.
            if ieVarInFile(variables, 'imgMean')
                disp('Saved using principal component method');
                load(fullname, 'imgMean')
                
                % Resample the image mean to the specified wavelength list
                if ~isempty(varargin)&& ~isempty(varargin{1})
                    extrapVal = 0;
                    imgMean = interp1(oldWave(:), imgMean(:), ...
                        newWave(:), 'linear', extrapVal); %#ok
                end
                
                % Sometimes we run out of memory here. So we should have a
                % try/catch sequence.
                %
                % The saved function was calculated using principal
                % components, not just the SVD. Hence, the mean is stored
                % and we must add it into the computed image.
                [photons, r, c] = RGB2XWFormat(photons);
                try
                    photons = repmat(imgMean(:), 1, r * c) + photons';
                catch ME
                    % Probably a memory error. Try with single precision.
                    if strcmp(ME.identifier, 'MATLAB:nomem')
                        photons = repmat(single(imgMean(:)), 1, r * c) ...
                            + single(photons');
                    else
                        ME.identifier
                    end
                end
                
                photons = double(XW2RGBFormat(photons', r, c));
                % figure(1);
                % imagesc(sum(img, 3));
                % axis image;
                % colormap(gray)
            else
                disp('Saved using svd method');
            end
            
            % Deal with the illuminant
            if ieVarInFile(variables, 'illuminant')
                load(fullname, 'illuminant')
            else
                % illuminant = [];
                warning('No illuminant information in %s\n', fullname);
            end
            
            % Force photons to be positive
            photons = max(photons, 0);
            
        else
            % The variable photons should be stored, there is no linear
            % model. We fill the basis slots. Also, we allow the photons
            % to be stored in 'photons' or 'data'. We allow the wavelength
            % to be stored in 'wave' or 'wavelength'.
            disp('Reading multispectral data with raw data.')
            
            
            if ieVarInFile(variables, 'photons')
                load(fullname, 'photons');
            elseif ieVarInFile(variables, 'data')
                load(fullname, 'data');
                photons = data;
                clear data;
            else
                error('No photon data in file');
            end
            if ieVarInFile(variables, 'comment')
                load(fullname, 'comment');
            end
            if ieVarInFile(variables, 'wave')
                load(fullname, 'wave');
            elseif ieVarInFile(variables, 'wavelength')
                load(fullname, 'wavelength');
                wave = wavelength;
                clear wavelength;
            end
            
            % Pull out the photons
            if ~isempty(varargin) && ~isempty(varargin{1})
                newWave = varargin{1};
                perfect = 0;
                idx = ieFindWaveIndex(wave, varargin{1}, perfect);
                photons = photons(:, :, idx);
                wave = newWave;
                % oldWave = wave;
                % wave = newWave;
            end
            basis.basis = [];
            basis.wave = round(wave);
        end
        
        % For linear model or no linear model, either way, we try to find
        % illuminant and resample.
        illuminant = [];
        if ieVarInFile(variables, 'illuminant')
            load(fullname, 'illuminant')
        else
            warning('No illuminant information in %s\n', fullname);
        end
        illuminant = illuminantModernize(illuminant);
        
        % Resample the illuminant to the specified wavelength list
        if ~isempty(varargin)&& ~isempty(varargin{1})
            % Resample the illuminant wavelength to the new wave in the
            % call to this function. This will also interpolate the
            % illuminant data.
            illuminant = illuminantSet(illuminant, 'wave', newWave(:));
        end
        
    otherwise
        fprintf('%s\n', imageType);
        error('Unknown image type.');
end

end