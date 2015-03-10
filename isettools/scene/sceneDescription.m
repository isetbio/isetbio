function txt = sceneDescription(scene)
%Text description of the scene properties, displayed in scene window
%
%    descriptionText = sceneDescription(scene)
%
% Copyright ImagEval Consultants, LLC, 2003.

if isempty(scene)
    txt = 'No scene';
else
    txt = sprintf('Name:\t%s\n',scene.name);
    r = sceneGet(scene,'rows'); c = sceneGet(scene,'cols');
    
    str = sprintf('(Row,Col):\t%.0f by %.0f \n',r,c);
    txt = addText(txt,str);
    
    u = round(log10(sceneGet(scene,'height','m')));
    if (u >= 0 ), units = 'm';
    elseif (u >= -3), units = 'mm'; 
    else units = 'um';
    end
    str = sprintf('Hgt,Wdth\t(%3.2f, %3.2f) %s\n', ...
                   sceneGet(scene,'height', units), ...
                   sceneGet(scene,'width', units), units);
    txt = addText(txt,str);

    u = round(log10(sceneGet(scene,'sampleSize', 'm')));
    if (u >= 0 ), units = 'm'; % 
    elseif (u >= -3), units = 'mm';
    else units = 'um';
    end
    str = sprintf('Sample:\t%3.2f %s \n', ...
                   sceneGet(scene,'sampleSize', units), units);
    txt = addText(txt, str);

    str = sprintf('Deg/samp: %2.2f\n',sceneGet(scene,'fov')/c);
    txt = addText(txt,str);
    
    wave = sceneGet(scene,'wave');
    spacing = sceneGet(scene,'binwidth');
    str = sprintf('Wave:\t%.0f:%.0f:%.0f nm\n', ...
                   min(wave(:)), spacing, max(wave(:)));
    txt = addText(txt, str);
    
    luminance = sceneGet(scene,'luminance');
    mx = max(luminance(:));
    mn = min(luminance(:));
    if mn == 0, 
        str = sprintf('DR: Inf\n  (max %.0f, min %.2f cd/m2)\n',mx,mn);
    else 
        dr = mx/mn;
        str = sprintf('DR: %.2f dB (max %.0f cd/m2)\n',20*log10(dr),mx);
    end
    
    txt = addText(txt,str);
    
end

end