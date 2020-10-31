function objNew = write(obj, varargin)
% Used by sceneEye.render. Typically not called directly.
%
% Syntax:
%   objNew = write(obj, [varargin])
%
% Description:
%	 A sceneEye object, has all the information needed to construct a PBRT
%	 file and render it. This function reads and interprets the parameters
%	 given in the sceneEye object and writes them out into an PBRT file.
%	 This file will later be rendered in sceneEye.render.
%
%    Copies the sceneEye parameters into the recipe and then writes out the
%    PBRT file. Typically a user will not run this directly, but rather it
%    will be run within the render function.
%
% Inputs:
%    obj   - The sceneEye object to render.  It contains the render recipe.
%
% Outputs:
%   objNew - Object. The object may have been modified in the processing
%            below. We return this modified version.
%
% Optional key/value pairs:
%    N/A
%
% See also
%   piWrite

%% PROGRAMMING TODO
%
%  This function should be using thisR.set/get not direct writes into the
%  recipe.
%
%  It seems to me also that this function copies the parameters in the
%  sceneEye object into the recipe object.  That is what a lot of the
%  sets/gets are about.  It then calls piWrite with the 'recipe'.  I am not
%  sure why sceneEye doesn't just manage the recipe directly, rather than
%  duplicating the parameters.
%

%% Make a copy of the current object
%
% We will render this copy, since we may make changes to certain parameters
% before rendering (i.e. in th eccentricity calculations) but we don't want
% these changes to show up original object given by the user.
objNew = copy(obj);
objNew.recipe = copy(obj.recipe);

%% Make some eccentricity calculations

% To render an image centered at a certain eccentricity without having
% change PBRT, we do the following:
% 1. Change the film size and resolution so that renders a larger image
%    that encompasses the desired eccentricity (tempWidth/tempHeight)
% 2. Insert a "crop window" PBRT parameter to only render the window
%    centered at the desired eccentricity with the desired film
%    diagonal/resolution.
ecc = objNew.eccentricity;

% I was having many bugs with my eccentricity code, so for now I've removed
% it for now. Ideally we do all the right calculations shown above and then
% use scene3d.recipe.set('cropwindow', [x1 x2 y1 y2]); and then carefully
% reset the angular support as well...
if(ecc ~= [0, 0])
    warning('Eccentricity is currently not implemented. Setting to zero.')
    ecc = [0 0];
end

%% Given the sceneEye object, make all adjustments needed to the recipe
thisR = objNew.recipe;

% Depending on the eye model, set the lens file appropriately
switch ieParamFormat(objNew.modelName)
    case {'navarro'}
        %{
        % Apply any accommodation changes
        if(isempty(objNew.accommodation))
            objNew.accommodation = 5;
            warning('No accommodation! Setting to 5 diopters.');
        end
        %}
        navarroWrite(objNew.recipe);
        
        % This function also writes out the Navarro lens file
        % recipe = setNavarroAccommodation(recipe, objNew.accommodation, ...
        %     objNew.workingDir);

    case {'legrand'}
        % Le Grand eye does not have accommodation (not yet at least).
        thisR = writeLegrandLensFile(thisR, objNew.workingDir); 

    case{'arizona'}
        if(isempty(objNew.accommodation))
            objNew.accommodation = 5;
            warning('No accommodation! Setting to 5 diopters.');
        end

        % This function also writes out the Arizona lens file.
        thisR = setArizonaAccommodation(thisR, objNew.accommodation, ...
            objNew.workingDir);

    case{'custom'}
        
        % Run this first to generate the IOR files.
        setNavarroAccommodation(thisR, 0, objNew.workingDir);

        % Copy the lens file given over
        if(isempty(obj.lensFile))
            error('No lens file given for custom eye.')
        else
            % Copy lens file over to the working directory and then attach
            % to recipe
            [success, message] = copyfile(obj.lensFile, objNew.workingDir);
            [~, n, e] = fileparts(obj.lensFile);

            if(success)
                thisR.camera.lensfile.value = ...
                    fullfile(objNew.workingDir, [n e]);
                thisR.camera.lensfile.type = 'string';
            else
                error('Error copying lens file. Err message: %s', message);
            end
        end
    case {'perspective'}
        % Probably in debug mode.  So let it go through, though you might
        % check debug mode.
    otherwise
        error('Unknown human eye model %s\n',modelName);

end

% semidiam = objNew.retinaDistance * tand(objNew.fov / 2);
if ~obj.debugMode
    % I am suspicious about the units here (mm vs m)
    semidiam = thisR.get('retina Distance') * tand(objNew.fov / 2);  % Looks like mm now
    thisR.set('retina semidiam',semidiam);  % Which makes this mm, too.
end

% recipe.camera.retinaSemiDiam.value = objNew.retinaDistance * ...
%    tand(objNew.fov / 2);

%{
if(strcmp(objNew.mmUnits, 'meters'))
    % We have some scenes that are in millimeters, not meters.  In that
    % case we set this flag to be true, indicating that mmUnits are true.
    thisR.camera.mmUnits.value = 'false';
    thisR.camera.mmUnits.type = 'bool';
end
%}

if(objNew.diffractionEnabled)
    % We should never get here.  Diffraction should always be set as
    % below.
    warning('Why am I here in diffraction enabled?');
    thisR.set('diffraction',true);
    % recipe.camera.diffractionEnabled.value = 'true';
    % recipe.camera.diffractionEnabled.type = 'bool';
end

% Renderer
numCABands = obj.recipe.get('num ca bands');
if(numCABands == 0 || numCABands == 1 || objNew.debugMode)
    % No spectral chromab rendering
    thisR.set('integrator subtype','path');
else
    numCABands = struct('value', objNew.numCABands, 'type', 'integer');
    thisR.set('integrator subtype','spectralpath');
    thisR.set('integrator num ca bands',numCABands);
end

%% Write out the adjusted recipe into a PBRT file
piWrite(thisR);

%{
pbrtFile = fullfile(objNew.workingDir, strcat(objNew.name, '.pbrt'));
recipe.set('outputFile', pbrtFile);
if(strcmp(recipe.exporter, 'C4D'))
    piWrite(recipe,  ...
        'overwritelensfile', false, 'overwriteresources', false, ...
        'creatematerials', true);
else
    piWrite(recipe,  ...
        'overwritelensfile', false, 'overwriteresources', false);
end
%}

% Not sure about any changes to recipe.  We should check carefully.
obj.recipe = thisR; 

end
