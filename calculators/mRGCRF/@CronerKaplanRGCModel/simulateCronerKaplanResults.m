function simulateCronerKaplanResults(obj, varargin)
    % Parse input
    p = inputParser;
    p.addParameter('synthesisOptions', obj.synthesisOptions, @isstruct);
    p.addParameter('generateVideo', true, @islogical);
    p.addParameter('generatePlots', true, @islogical);
    p.parse(varargin{:});
    generateVideo = p.Results.generateVideo;
    generatePlots = p.Results.generatePlots;
    synthesisOptions = p.Results.synthesisOptions;
    
    % ecc range appropriate for the Croner & Kaplan '95 study   
    examinedEccDegs = abs(normrnd(5,10,[1 numel(obj.centerData('size').eccDegs)]));
    
    
    if (generateVideo)
        videoOBJ1 = VideoWriter('params', 'MPEG-4'); % H264 format
        videoOBJ1.FrameRate = 2;
        videoOBJ1.Quality = 100;
        videoOBJ1.open();

        videoOBJ2 = VideoWriter('ratios', 'MPEG-4'); % H264 format
        videoOBJ2.FrameRate = 2;
        videoOBJ2.Quality = 100;
        videoOBJ2.open();

        videoOBJ3 = VideoWriter('rf2d', 'MPEG-4'); % H264 format
        videoOBJ3.FrameRate = 2;
        videoOBJ3.Quality = 100;
        videoOBJ3.open();
        
        videoOBJ4 = VideoWriter('rf1d', 'MPEG-4'); % H264 format
        videoOBJ4.FrameRate = 2;
        videoOBJ4.Quality = 100;
        videoOBJ4.open();
        
        trialsNum = 100;
    else
        trialsNum = 1;
    end
    
    for trialNo = 1:trialsNum   
        obj.synthesizeData(examinedEccDegs, synthesisOptions);
        [hFig1, hFig2, hFig3, hFig4] = obj.plotSynthesizedData();
        drawnow;
        if (generateVideo)
            videoOBJ1.writeVideo(getframe(hFig1));
            videoOBJ2.writeVideo(getframe(hFig2));
            videoOBJ3.writeVideo(getframe(hFig3));
            videoOBJ4.writeVideo(getframe(hFig4));
        end
    end
    
    if (generateVideo)
        videoOBJ1.close();
        videoOBJ2.close();
        videoOBJ3.close();
        videoOBJ4.close();
    end
    
    if (generatePlots)
        obj.plotlabOBJ.exportFig(hFig1, 'pdf', 'synthesizedParams', pwd());
        obj.plotlabOBJ.exportFig(hFig2, 'pdf', 'synthesizedRatios', pwd());
        obj.plotlabOBJ.exportFig(hFig3, 'pdf', 'synthesized2DRFs', pwd());
    end
end

