function thisR = setNavarroAccommodation(thisR, accommodation, workingFolder)
% Deprecated:  
% 
% Change lens model to match the accommodation
%
% Syntax:
%   thisR = setNavarroAccommodation(thisR, accommodation, [workingFolder])
%
% Description:
%    For the human eye models, when accommodation changes the lens
%    file changes. We write these new files out and reference them in
%    the structure. 
%
%   This scope of the Navarro model includes accommodation from 0 to
%   10 diopters. A value of 0 diopters means the focus is at infinity.
%   10 diopters means the eye model focus is at 0.1 meter. 
%
%   However, the Navarro accommodation values do not match the Zemax
%   values.  So we call a function that converts the user's
%   accommodation value to the one from Zemax.  We should probably
%   allow ourselves to turn this conversion off (BW).
%
% Inputs:
%    thisR         - Object (Render recipe)
%    accommodation - Numeric. The accommodation to shape the modified
%                    thisR by.
%
% Optional
%    workingFolder - String. By default it is
%                    thisR.get('lens output dir')
%
% Outputs:
%    thisR          - The modified render recipe.
%
%
% See also
%   navarroLensCreate

% History:
%    xx/xx/17  TL   Created by Trisha Lian IESTBIO Team 2017
%    12/19/17  jnm  Formatting
%    01/25/19  JNM  Changed output of filenames to remove extra periods,
%                   which is breaking Windows executions.
%    05/29/19  JNM  Second documentation pass (minor tweaks)

%% Check and make sure this recipe has a human eye model
if ~strcmp(thisR.get('camera subtype'),'humaneye')
    warning('The camera type is not a human eye model. Returning untouched.');
    return;
end

%% Check inputs
if(~(accommodation >= 0 && accommodation <= 10))
    % This is the scope of the model.  0 diopters means the focus is at
    % infinity.  10 diopters means the eye model focus is at 0.1 meter.
    error('Accommodation must be between 0 and 10 diopters.');
end

% Simply over-write the eye model with the new accommodation
navarroWrite(thisR,accommodation);

end



