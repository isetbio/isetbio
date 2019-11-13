function generateAndDockAllFigures(obj)

    desktop = com.mathworks.mde.desk.MLDesktop.getInstance;
    WatsonRGCpaperFiguresGroup = desktop.addGroup('WatsonRGCpaperFiguresGroup');
    desktop.setGroupDocked('WatsonRGCpaperFiguresGroup', 0);
    groupedFigsCols = 3;
    groupedFigsRows = 3;
    GroupDims = java.awt.Dimension(groupedFigsCols, groupedFigsRows);
    % 1: Maximized, 2: Tiled, 3: Floating
    desktop.setDocumentArrangement('WatsonRGCpaperFiguresGroup', 2, GroupDims)
    hFigs = gobjects(1, groupedFigsCols*groupedFigsRows);
    bakWarn = warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
    
    for iFig = 1:groupedFigsCols*groupedFigsRows
        hFigs(iFig) = figure('WindowStyle', 'docked', ...
                'Name', sprintf('Figure %d', iFig), ...
                'NumberTitle', 'off');
        
        pause(0.1);  % Magic, reduces rendering errors
        set(get(handle(hFigs(iFig)), 'javaframe'), 'GroupName', 'WatsonRGCpaperFiguresGroup');
    end
    
    for iFig = 1:groupedFigsCols*groupedFigsRows
        switch(iFig)
            case 1
                % Cone density as a function of eccentricity for all quadrants
                obj.generateFigure1(hFigs(iFig));
            case 2
                % RF density of all RGCs as a function of eccentricity for all quadrants
                obj.generateFigure5(hFigs(iFig));
            case 3
                % fraction of midget to total RGCs RFs as a function of eccentrity
                obj.generateFigure8(hFigs(iFig));
            case 4
                % RF density of midget RGCs as a function of eccentricity for all quadrants
                obj.generateFigure9(hFigs(iFig));
            case 5
                % fraction of midget to total RGCs RFs as a function of eccentrity
                obj.generateFigure10(hFigs(iFig));
            case 6
                % fraction of midget to total RGCs RFs as a function of eccentrity
                obj.generateFigure11(hFigs(iFig));
            case 7
                % Ratio of midget RGCs to cones as a function of eccentricity for all quadrants
                obj.generateFigure14(hFigs(iFig));
            case 8
                % Relation between retinal distance in mm and degs
                obj.generateFigureA1(hFigs(iFig));
        end
        drawnow;
        pause(0.1);
        hFigs(iFig);
    end
    warning(bakWarn);

end

