function setWave(obj, src)
% Keep photopigment and macular pigment wavelength sampling consistent.
%
% Syntax:
%   setWave(obj, src)
%
% Description:
%    When you set the photpigment wavelength sampling, it sets the macular
%    pigment wavelength sampling and vice versa.
%
% Inputs:
%    obj - The cone mosaic object
%    src - The light source? 
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note: DHB - THIS MAKES ME NERVOUS.  NEED A FULLER EXPLANATION AND
%      TO CHECK THROUGH THE LOGIC. THIS RELATES TO THE BROAD QUESTION OF
%      WHO IS IN CHARGE OF KEEPING WAVELENGTH SAMPLING CONSISTENT.]
%    * [Note: DHB - I DON'T UNDERSTAND THE PURPOSE OF THE ~ AS THE LAST
%      INPUT ARGUMENT.]
%

% History:
%    xx/xx/16  HJ   ISETBIO Team 2016
%    02/23/18  jnm  Formatting

switch src.DefiningClass.Name
    case 'photoPigment'
        obj.macular.wave = obj.pigment.wave;
    case 'Macular'
        obj.pigment.wave = obj.macular.wave;
end
end