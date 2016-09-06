function setWave(obj, src, ~)
% callback function for listeners
% When you set the cone, it sets the macular, and vice versa.

switch src.DefiningClass.Name
    case 'photoPigment'
        obj.macular.wave = obj.pigment.wave;
    case 'Macular'
        obj.pigment.wave = obj.macular.wave;
end
end