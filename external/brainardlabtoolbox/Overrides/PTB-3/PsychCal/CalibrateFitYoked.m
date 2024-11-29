function [varargout] = CalibrateFitYoked(calOrCalStruct)
% Usage (new style):
% CalibrateFitYoked(calOrCalStruct)
% where calOrCalStruct is a @CalStruct object.
%
% Usage (old style):
% cal = CalibrateFitYoked(calOrCalStruct)
% where cal and calOrCalStruct are both old-style cal structures.
%
% Fit the gamma data from the yoked measurements.
%
% This has to do with Brainard lab HDR display calibration procedures
% and doesn't do anything unless some special fields exist
% in the calibration structure.  It's in the PTB because when we use
% it, we want to call it from script RefitCalGamma, and that one does
% belong in the PTB.  
%
% 4/30/10  dhb, kmo, ar Wrote it.
% 5/24/10  dhb          Update comment.
% 5/25/10  dhb, ar      New input format.
% 5/28/10  dhb, ar      Execute conditionally.
% 6/10/10  dhb          Make sure returned gamma values in range an monotonic.
% 5/26/12  dhb          Improve the comment so this is a little less weird.
% 5/28/14  npc          Modifications for accessing calibration data using a @CalStruct object.
%                       The first input argument can be either a @CalStruct object (new style), or a cal structure (old style).
%                       Passing a @CalStruct object is the preferred way because it results in 
%                       (a) less overhead (@CalStruct objects are passed by reference, not by value), and
%                       (b) better control over how the calibration data are accessed.

    %% Debugging switch
    DEBUG = 0;

    % Specify @CalStruct object that will handle all access to the calibration data.
    [calStructOBJ, inputArgIsACalStructOBJ] = ObjectToHandleCalOrCalStruct(calOrCalStruct);
    
    if (inputArgIsACalStructOBJ)
        % The input (calOrCalStruct) is a @CalStruct object. Make sure that CalibrateFitYoked is called with no return variables.
        if (nargout > 0)
            error('There should be NO return parameters when calling CalibrateFitYoked with a @CalStruct input. For more info: doc CalibrateFitYoked.');
        end
    else
        % The input (calOrCalStruct) is a cal struct. Clear it to avoid  confusion.
        clear 'calOrCalStruct';
    end
    % From this point onward, all access to the calibration data is accomplised via the calStructOBJ.


    % Extract needed data
    yokedMethod     = calStructOBJ.get('yokedmethod');
    yokedSpectra    = calStructOBJ.get('yoked.spectra');
    yokedSettings   = calStructOBJ.get('yoked.settings');
    nPrimaryBases   = calStructOBJ.get('nPrimaryBases');
    P_device        = calStructOBJ.get('P_device');

    %% Check that this is possible
    OKTODO = 1;
    if ((isempty(yokedMethod)) || ((~isempty(yokedMethod))&&(yokedMethod==0))) 
        OKTODO = 0;
    end
    if (isempty(yokedSpectra))
        OKTODO = 0;
    end
    if (nPrimaryBases == 0)
        OKTODO = 0;
    end
    if (~OKTODO)
        if (~inputArgIsACalStructOBJ)
            % Old-style functionality. Return modified cal.
            varargout{1} = calStructOBJ.cal;
        end
        return;
    end

    %% Average yoked measurements for this primary
    yokedSpds = yokedSpectra;

    %% Fit each spectrum with the linear model for all three primaries
    % and project down onto this
    projectedYokedSpd = P_device*(P_device\yokedSpds);

    %% Now we have to adjust the linear model so that it has our standard 
    % properties.

    % Make first three basis functions fit maxSpd exactly
    maxSpd = projectedYokedSpd(:,end);
    weights = P_device\maxSpd;
    currentLinMod = zeros(size(P_device));
    for i = 1:nDevices
        tempLinMod = 0;
        for j = 1:nPrimaryBases
            tempLinMod = tempLinMod + P_device(:,i+(j-1)*nDevices)*weights(i+(j-1)*nDevices);
        end
        currentLinMod(:,i) = tempLinMod;
    end
    weights = currentLinMod(:,1:nDevices)\maxSpd;
    for i = 1:nDevices
        currentLinMod(:,i) = currentLinMod(:,i)*weights(i);
    end
    maxPow = max(max(currentLinMod(:,1:nDevices)));

    % Now find the rest of the linear model
    clear tempLinMod
    for i = 1:nDevices
        for j = 1:nPrimaryBases
            tempLinMod(:,j) = P_device(:,i+(j-1)*nDevices); %#ok<AGROW>
        end
        residual = tempLinMod - currentLinMod(:,i)*(currentLinMod(:,i)\tempLinMod);
        restOfLinMod = FindLinMod(residual,nPrimaryBases-1);
        for j = 2:nPrimaryBases
            tempMax = max(abs(restOfLinMod(:,j-1)));
            currentLinMod(:,i+(j-1)*nDevices) = maxPow*restOfLinMod(:,j-1)/tempMax;
        end
    end

    % Span of P_device and currentLinMod should be the same.  Check this.
    if (DEBUG)
        check = currentLinMod - P_device*(P_device\currentLinMod);
        if (max(abs(check(:))) > 1e-10)
            error('Two linear models that should have the same span don''t');
        end
    end

    % Express yoked spectra in terms of model weights
    gammaTable = currentLinMod\yokedSpectra;
    tempSpd = currentLinMod*gammaTable;
    for i = 1:nDevices
        index = gammaTable(i,:) > 1;
        gammaTable(i,index) = 1;
        gammaTable(i,:) = MakeMonotonic(HalfRect(gammaTable(i,:)'))';
    end

    % Stash info in calibration structure
    calStructOBJ.set('P_device', currentLinMod);

    % When R=G=B, we just use the common settings.
    if (yokedmethod == 1)
        calStructOBJ.set('rawGammaInput', yokedSettings(1,:)');
        calStructOBJ.set('rawGammaTable', gammaTable');
    % When measurements are at a specified chromaticity, need to interpolate gamma
    % functions so that we have them for each device on a common scale.
    elseif (yokedmethod == 2)
        calStructOBJ.set('rawGammaInput', yokedSettings');
        calStructOBJ.set('rawGammaTable', gammaTable');
    end

    %% Debugging
    if (DEBUG)
        S = calStructOBJ.get('S');
        load T_xyz1931
        T_xyz=683*SplineCmf(S_xyz1931, T_xyz1931, S);

        % Meausured xyY
        measuredYokedxyY = XYZToxyY(T_xyz*yokedSpectra);

        % Raw linear model fit xyY
        projectedYokedxyY = XYZToxyY(T_xyz*projectedYokedSpd);

        % Predicted xyY
        predictedSpd = P_device*rawGammaTable';
        predictedYokedxyY = XYZToxyY(T_xyz*predictedSpd);

        % Plot luminance obtained vs. desired
        [lumPlot,f] = StartFigure('standard');
        f.xrange = [0 size(yokedSettings, 2)]; f.nxticks = 6;
        f.yrange = [0 360]; f.nyticks = 5;
        f.xtickformat = '%0.0f'; f.ytickformat = '%0.2f ';
        plot(measuredYokedxyY(3,:)','ro','MarkerSize',f.basicmarkersize,'MarkerFaceColor','r');
        plot(projectedYokedxyY(3,:)','bo','MarkerSize',f.basicmarkersize,'MarkerFaceColor','b');
        plot(predictedYokedxyY(3,:)','go','MarkerSize',f.basicmarkersize,'MarkerFaceColor','g');

        xlabel('Test #','FontName',f.fontname,'FontSize',f.labelfontsize);
        ylabel('Luminance (cd/m2)','FontName',f.fontname,'FontSize',f.labelfontsize);
        FinishFigure(lumPlot,f);

        % Plot x chromaticity obtained vs. desired
        [xPlot,f] = StartFigure('standard');
        f.xrange = [0 size(yokedSettings, 2)]; f.nxticks = 6;
        f.yrange = [0.2 0.6]; f.nyticks = 5;
        f.xtickformat = '%0.0f'; f.ytickformat = '%0.2f ';
        plot(measuredYokedxyY(1,:)','ro','MarkerSize',f.basicmarkersize,'MarkerFaceColor','r');
        plot(projectedYokedxyY(1,:)','bo','MarkerSize',f.basicmarkersize,'MarkerFaceColor','b');
        plot(predictedYokedxyY(1,:)','go','MarkerSize',f.basicmarkersize,'MarkerFaceColor','g');
        xlabel('Test #','FontName',f.fontname,'FontSize',f.labelfontsize);
        ylabel('x chromaticity','FontName',f.fontname,'FontSize',f.labelfontsize);
        FinishFigure(xPlot,f);

        % Plot y chromaticity obtained vs. desired
        [yPlot,f] = StartFigure('standard');
        f.xrange = [0 size(yokedSettings, 2)]; f.nxticks = 6;
        f.yrange = [0.2 0.6]; f.nyticks = 5;
        f.xtickformat = '%0.0f'; f.ytickformat = '%0.2f ';
        plot(measuredYokedxyY(2,:)','ro','MarkerSize',f.basicmarkersize,'MarkerFaceColor','r');
        plot(projectedYokedxyY(2,:)','bo','MarkerSize',f.basicmarkersize,'MarkerFaceColor','b');
        plot(predictedYokedxyY(2,:)','go','MarkerSize',f.basicmarkersize,'MarkerFaceColor','g');
        xlabel('Test #','FontName',f.fontname,'FontSize',f.labelfontsize);
        ylabel('y chromaticity','FontName',f.fontname,'FontSize',f.labelfontsize);
        FinishFigure(yPlot,f);

        drawnow;
    end

    if (~inputArgIsACalStructOBJ)
        % Old-style functionality. Return modified cal.
        varargout{1} = calStructOBJ.cal;
        % calStructOBJ is not needed anymore. So clear it from the memory.
        clear 'calStructOBJ'
    end

end



