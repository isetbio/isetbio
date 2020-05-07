function displayPSFs()

    goodSubjects = [8 9 10];
    
    % Compute good subject PSFs along the Y = 0, negative X, with an
    % eccentricity resolution of 1 deg
    eccXrange = [-40 0];
    eccYrange = [0 0];
    deltaEcc = 1;
    
    % Compute the PSFs at the desired eccentricities
    [hEcc, vEcc, thePSFs, thePSFsupportDegs] = CronerKaplanRGCModel.psfAtEccentricity(goodSubjects, ...
        eccXrange, eccYrange, deltaEcc);
    
    % Visualize them
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 2040 300]);

    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', numel(goodSubjects), ...
           'colsNum', numel(hEcc), ...
           'heightMargin',  0.001, ...
           'widthMargin',    0.005, ...
           'leftMargin',     0.02, ...
           'rightMargin',    0.00, ...
           'bottomMargin',   0.02, ...
           'topMargin',      0.01);
       
    cMap = brewermap(512, 'greys');
    for subjIdx = 1:numel(goodSubjects)
       for eccYIndex = 1:numel(vEcc)
       for eccXIndex = 1:numel(hEcc)
            
            ax = subplot('Position', subplotPosVectors(subjIdx, eccXIndex).v);
            imagesc(ax, thePSFsupportDegs, thePSFsupportDegs, squeeze(thePSFs(subjIdx, eccYIndex, eccXIndex,:,:)));
            axis(ax, 'square'); 
            set(ax, 'XTickLabel', {}, 'YTickLabel', {});
            if ((vEcc(eccYIndex) == 0) && (hEcc(eccXIndex) == 0))
                title(ax,sprintf('%2.0f at %2.0f^o', goodSubjects(subjIdx), hEcc(eccXIndex)));
            end
            
            colormap(ax,cMap);
            drawnow;
       end
       end
    end % subIdx
    
    
    % Custom 2D sampling (2.5 deg) for one subject
    eccXrange = [-20 20];
    eccYrange = [-20 20];
    deltaEcc = 2.5;
    goodSubjects = [8];
    
    % Compute the PSFs at the desired eccentricities
    [hEcc, vEcc, thePSFs, thePSFsupportDegs] = CronerKaplanRGCModel.psfAtEccentricity(goodSubjects, ...
        eccXrange, eccYrange, deltaEcc);
    
    hFig = figure(2); clf;
    set(hFig, 'Position', [10 10 1200 1180]);

    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', numel(vEcc), ...
           'colsNum', numel(hEcc), ...
           'heightMargin',  0.01, ...
           'widthMargin',   0.001, ...
           'leftMargin',    0.02, ...
           'rightMargin',   0.01, ...
           'bottomMargin',  0.02, ...
           'topMargin',     0.01);

     for eccYIndex = 1:numel(vEcc)
     for eccXIndex = 1:numel(hEcc)
        ax = subplot('Position', subplotPosVectors(eccYIndex, eccXIndex).v);
        imagesc(ax, thePSFsupportDegs, thePSFsupportDegs, squeeze(thePSFs(1, eccYIndex, eccXIndex,:,:)));
        axis(ax, 'square'); 
        set(ax, 'XTickLabel', {}, 'YTickLabel', {});
        title(ax,sprintf('%2.1f,%2.1f', hEcc(eccXIndex), vEcc(eccYIndex)));
        colormap(ax,cMap);
        drawnow;
     end
     end
     
end

