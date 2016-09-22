function t_coneMosaicDemosaicResponses

    % Set up scene
    scene = sceneCreate; 
    scene = sceneSet(scene, 'h fov', 1.0);
    oi = oiCreate('human'); 
    oi = oiCompute(oi, scene);
   
    % Set up mosaic
    cMosaicOBJ = coneMosaic();
    cMosaicOBJ.setSizeToFOV([sceneGet(scene, 'h fov'), sceneGet(scene, 'v fov')]);
    cMosaicOBJ.noiseFlag = false;
    [~, currentsMap] = cMosaicOBJ.compute(oi,'currentFlag', true);

    % Call demosaicing methods
    [demosaicedIsomerizationsMaps, sRGB] = cMosaicOBJ.demosaicedIsomerizationMaps();
    demosaicedPhotoCurrentMaps = cMosaicOBJ.demosaicedPhotoCurrentMaps(currentsMap);
  
    % Display results
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 1340 470]);
    cLims = [min(demosaicedIsomerizationsMaps(:)) max(demosaicedIsomerizationsMaps(:))];
    subplot(2,3,1);
    imagesc(squeeze(demosaicedIsomerizationsMaps(:,:,1,1)), cLims);
    title('Lcone isomerization demosaiced map');
    colorbar()
    axis 'image';
   
    subplot(2,3,2);
    imagesc(squeeze(demosaicedIsomerizationsMaps(:,:,2,1)), cLims);
    title('Mcone isomerization demosaiced map');
    colorbar()
    axis 'image';
   
    subplot(2,3,3);
    imagesc(squeeze(demosaicedIsomerizationsMaps(:,:,3,1)), cLims);
    title('Scone isomerization demosaiced map');
    colorbar()
    axis 'image';
  

    cLims = [min(demosaicedPhotoCurrentMaps(:)) max(demosaicedPhotoCurrentMaps(:))];
    subplot(2,3,4);
    imagesc(squeeze(demosaicedPhotoCurrentMaps(:,:,1,1)), cLims);
    title('Lcone photocurrent demosaiced map');
    colorbar()
    axis 'image';
   
    subplot(2,3,5);
    imagesc(squeeze(demosaicedPhotoCurrentMaps(:,:,2,1)), cLims);
    title('Mcone photocurrent demosaiced map');
    colorbar()
    axis 'image';
   
    subplot(2,3,6);
    imagesc(squeeze(demosaicedPhotoCurrentMaps(:,:,3,1)), cLims);
    title('Scone photocurrent demosaiced map');
    colorbar()
    axis 'image';
   
    colormap(gray);
    drawnow;
   
    hFig = figure(2); clf;
    set(hFig, 'Position', [10 10 1470 990])
    subplot(2,2,1);
    imshow(sceneGet(scene, 'RGB'));
    title('input scene RGB rendition');
   
    subplot(2,2,2);
    imshow(sRGB(:,:,:,1));
    title('photoisomerization RGB rendition');
    
    subplot(2,2,4);
    uData = cMosaicOBJ.plot('cone mosaic', 'hf', 'none');
    imagesc(uData.mosaicImage); axis off; axis image;
    title('cone mosaic');
end

