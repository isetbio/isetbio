function setWave(obj, src, ~)
% When you set the cone, it sets the macular, and vice versa.
% 
% HJ ISETBIO Team 2016

switch src.DefiningClass.Name
    case 'photoPigment'
        obj.macular.wave = obj.pigment.wave;
    case 'Macular'
        obj.pigment.wave = obj.macular.wave;
end
end