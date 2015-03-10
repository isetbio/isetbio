function plotAbso(absorptions)

warning('Deprecated. Use sensorVisualize instead');
if notDefined('absorptions.data')
    error('No data field in absorptions')
end

frameNumber = floor(rand*size(absorptions.data,3))+1;

if isfield(absorptions, 'unadapted')
    figure; subplot(121);
    imagesc(absorptions.unadapted(:,:,frameNumber), [0, 2*mean(mean(absorptions.unadapted(:,:,frameNumber)))]); 
    colormap(gray)
    colorbar
    title('Absorptions before adaptation')

    subplot(122)
    imagesc(absorptions.data(:,:,frameNumber), [0, 2*mean(mean(absorptions.data(:,:,frameNumber)))]); colormap(gray)
    if isfield(absorptions, 'typeAdapt')
        title(sprintf('Absorptions after adaptation (type %d)', absorptions.typeAdapt))
    else
        title('Absorptions after adaptation')   
    end
    colorbar
    
else
    figure;
    imagesc(absorptions.data(:,:,frameNumber), [0, 2*mean(mean(absorptions.data(:,:,frameNumber)))]); colormap(gray)
    title('Absorptions')
    colorbar
end

end
