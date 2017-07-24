function setWave(obj, src, ~)
%SETWAVE  Keep photopigment and macular pigment wavelength sampling consistent.
%
%   When you set the photpigment wavelength sampling, it sets the macular pigment wavelength sampling
%   and vice versa.
%
%   [DHB COMMENT: THIS MAKES ME NERVOUS.  NEED A FULLER EXPLANATION AND TO CHECK THROUGH THE LOGIC.  
%   THIS RELATES TO THE BROAD QUESTION OF WHO IS IN CHARGE OF KEEPING WAVELENGTH SAMPLING CONSISTENT.]
%   [DHB COMMENT: I DON'T UNDERSTAND THE PURPOSE OF THE ~ AS THE LAST INPUT ARGUMENT.]

% HJ ISETBIO Team 2016

switch src.DefiningClass.Name
    case 'photoPigment'
        obj.macular.wave = obj.pigment.wave;
    case 'Macular'
        obj.pigment.wave = obj.macular.wave;
end
end