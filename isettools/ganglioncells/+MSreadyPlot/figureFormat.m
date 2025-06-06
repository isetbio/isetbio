function f = figureFormat(panelLayout)

    f.fontSize = 28;
    f.markerSize = 16;
    f.lineWidth = 2.0;
    f.LconeColor = [1 0.2 0.6];
    f.MconeColor = [0.2 1 0.4];
    f.SconeColor = [0.6 0.2 0.9];
    f.axisLineWidth = 1.0;

    f.axisColor = [0.3 0.3 0.3];
    f.axisFontAngle = 'italic';

    f.titleFontSize = 20;
    f.titleColor = [0.2 0.2 0.2];
    f.titleFontWeight = 'normal';

    f.legendFontSize = 24;

    f.axisOffsetFactor = -0.03;

    % figure sizes
    switch panelLayout
        case '1x4'
            f.figureSize = [1500 450];
            f.subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                'rowsNum', 1, ...
                'colsNum', 4, ...
                'heightMargin',  0.01, ...
                'widthMargin',    0.06, ...
                'leftMargin',     0.04, ...
                'rightMargin',    0.01, ...
                'bottomMargin',   0.15, ...
                'topMargin',      0.05);
        case '2x4'
            f.figureSize = [1500 850];
            f.subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                'rowsNum', 2, ...
                'colsNum', 4, ...
                'heightMargin',  0.08, ...
                'widthMargin',    0.06, ...
                'leftMargin',     0.05, ...
                'rightMargin',    0.00, ...
                'bottomMargin',   0.06, ...
                'topMargin',      0.02);

        case '2x4-tall'
            f.figureSize = [2000 1160];
            f.subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                'rowsNum', 2, ...
                'colsNum', 4, ...
                'heightMargin',  0.08, ...
                'widthMargin',    0.06, ...
                'leftMargin',     0.05, ...
                'rightMargin',    0.00, ...
                'bottomMargin',   0.06, ...
                'topMargin',      0.02);
            
        case '4x4'
            f.figureSize = [1200 1150];
            f.subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                'rowsNum', 4, ...
                'colsNum', 4, ...
                'heightMargin',  0.10, ...
                'widthMargin',    0.06, ...
                'leftMargin',     0.05, ...
                'rightMargin',    0.00, ...
                'bottomMargin',   0.06, ...
                'topMargin',      0.02);

       case '2x3'
            f.figureSize = [2025 1200];
            f.subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                'rowsNum', 2, ...
                'colsNum', 3, ...
                'heightMargin',  0.04, ...
                'widthMargin',    0.04, ...
                'leftMargin',     0.03, ...
                'rightMargin',    0.00, ...
                'bottomMargin',   0.02, ...
                'topMargin',      0.01);

       case '3x3'
            f.figureSize = [950 975]*1.25;
            f.subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                'rowsNum', 3, ...
                'colsNum', 3, ...
                'heightMargin',  0.05, ...
                'widthMargin',    0.04, ...
                'leftMargin',     0.08, ...
                'rightMargin',    0.00, ...
                'bottomMargin',   0.07, ...
                'topMargin',      0.01);

       case '4x8'
            f.figureSize = [2050 1150];
            f.subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                'rowsNum', 4, ...
                'colsNum', 8, ...
                'heightMargin',  0.04, ...
                'widthMargin',    0.03, ...
                'leftMargin',     0.04, ...
                'rightMargin',    0.00, ...
                'bottomMargin',   0.05, ...
                'topMargin',      0.03);

        case '2x2'
            f.figureSize = [800 800];
            f.subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                'rowsNum', 2, ...
                'colsNum', 2, ...
                'heightMargin',  0.12, ...
                'widthMargin',    0.08, ...
                'leftMargin',     0.08, ...
                'rightMargin',    0.01, ...
                'bottomMargin',   0.07, ...
                'topMargin',      0.02);

        case '1x3 RF poster'
            f.figureSize = [1080 390];
            f.fontSize = 20;
            f.subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                'rowsNum', 1, ...
                'colsNum', 3, ...
                'heightMargin',  0.12, ...
                'widthMargin',    0.06, ...
                'leftMargin',     0.07, ...
                'rightMargin',    0.01, ...
                'bottomMargin',   0.10, ...
                'topMargin',      0.00);

        case '1x4 RF poster'
            f.figureSize = [1300 350];
            f.fontSize = 20;
            f.markerSize = f.markerSize - 6;
            f.legendFontSize = f.legendFontSize - 6;
            f.subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                'rowsNum', 1, ...
                'colsNum', 4, ...
                'heightMargin',  0.12, ...
                'widthMargin',    0.055, ...
                'leftMargin',     0.065, ...
                'rightMargin',    0.01, ...
                'bottomMargin',   0.09, ...
                'topMargin',      0.00);


        case '1x1 small'
            % Used in PLOS One 2023 paper figures
            f.figureSize = [700 500];
            f.markerSize = 20;
            f.legendFontSize = 20;
            f.subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                'rowsNum', 1, ...
                'colsNum', 1, ...
                'heightMargin',  0.00, ...
                'widthMargin',    0.00, ...
                'leftMargin',     0.15, ...
                'rightMargin',    0.03, ...
                'bottomMargin',   0.18, ...
                'topMargin',      0.05);

        case '1x1 small no labels'
            % Used in PLOS One 2023 paper figures
            f.figureSize = [500 500];
            f.markerSize = 20;
            f.legendFontSize = 20;
            f.subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                'rowsNum', 1, ...
                'colsNum', 1, ...
                'heightMargin',  0.00, ...
                'widthMargin',    0.00, ...
                'leftMargin',     0.01, ...
                'rightMargin',    0.01, ...
                'bottomMargin',   0.01, ...
                'topMargin',      0.01);

        case '1x1 medium'
            % Used in PLOS One 2023 paper figures
            f.figureSize = [1280 1280];
            f.markerSize = 20;
            f.legendFontSize = 20;
            f.subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                'rowsNum', 1, ...
                'colsNum', 1, ...
                'heightMargin',  0.00, ...
                'widthMargin',    0.00, ...
                'leftMargin',     0.06, ...
                'rightMargin',    0.03, ...
                'bottomMargin',   0.06, ...
                'topMargin',      0.0);

        case '1x1 small tall'
            f.figureSize = [700 900];
            f.subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                'rowsNum', 1, ...
                'colsNum', 1, ...
                'heightMargin',  0.00, ...
                'widthMargin',    0.00, ...
                'leftMargin',     0.15, ...
                'rightMargin',    0.03, ...
                'bottomMargin',   0.097, ...
                'topMargin',      0.025);

        case '1x1 small wide'
            f.figureSize = [1600 700];
            f.subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                'rowsNum', 1, ...
                'colsNum', 1, ...
                'heightMargin',  0.00, ...
                'widthMargin',    0.00, ...
                'leftMargin',     0.04, ...
                'rightMargin',    0.01, ...
                'bottomMargin',   0.10, ...
                'topMargin',      0.01);

        case '1x1 small ultra wide'
            f.figureSize = [2000 750];
            f.subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                'rowsNum', 1, ...
                'colsNum', 1, ...
                'heightMargin',  0.00, ...
                'widthMargin',    0.00, ...
                'leftMargin',     0.05, ...
                'rightMargin',    0.01, ...
                'bottomMargin',   0.12, ...
                'topMargin',      0.02);

        case '1x1 small wide tall'
            f.figureSize = [1500 1300];
            f.subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                'rowsNum', 1, ...
                'colsNum', 1, ...
                'heightMargin',  0.00, ...
                'widthMargin',    0.00, ...
                'leftMargin',     0.10, ...
                'rightMargin',    0.02, ...
                'bottomMargin',   0.13, ...
                'topMargin',      0.00);

        case '1x1 small very wide'
            f.figureSize = [1000 500];
            f.subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                'rowsNum', 1, ...
                'colsNum', 1, ...
                'heightMargin',  0.00, ...
                'widthMargin',    0.00, ...
                'leftMargin',     0.06, ...
                'rightMargin',    0.01, ...
                'bottomMargin',   0.14, ...
                'topMargin',      0.00);

        case '1x1 large'
            f.figureSize = [800 800];
            f.subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                'rowsNum', 1, ...
                'colsNum', 1, ...
                'heightMargin',  0.00, ...
                'widthMargin',    0.00, ...
                'leftMargin',     0.08, ...
                'rightMargin',    0.00, ...
                'bottomMargin',   0.07, ...
                'topMargin',      0.02);

       case '1x1 long'
            f.figureSize = [1200 650];
            f.subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                'rowsNum', 1, ...
                'colsNum', 1, ...
                'heightMargin',  0.00, ...
                'widthMargin',    0.00, ...
                'leftMargin',     0.07, ...
                'rightMargin',    0.01, ...
                'bottomMargin',   0.10, ...
                'topMargin',      0.02);

        case '1x1 very long poster'
            f.fontSize = 24;
            f.figureSize = [2050 250];
            f.subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                'rowsNum', 1, ...
                'colsNum', 1, ...
                'heightMargin',  0.00, ...
                'widthMargin',    0.00, ...
                'leftMargin',     0.04, ...
                'rightMargin',    0.00, ...
                'bottomMargin',   0.09, ...
                'topMargin',      -0.03);
            f.tickLength = [0.01/4 0.01/10];

        case '1x1 poster'
            f.fontSize = 24;
            f.figureSize = [2050 1150];
            f.subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                'rowsNum', 1, ...
                'colsNum', 1, ...
                'heightMargin',  0.00, ...
                'widthMargin',    0.00, ...
                'leftMargin',     0.05, ...
                'rightMargin',    0.00, ...
                'bottomMargin',   0.01, ...
                'topMargin',      -0.06);
            f.tickLength = [0.01/4 0.01/10];

        case '1x2 large'
            f.figureSize = [1700 800];
            f.subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                'rowsNum', 1, ...
                'colsNum', 2, ...
                'heightMargin',  0.00, ...
                'widthMargin',    0.06, ...
                'leftMargin',     0.08, ...
                'rightMargin',    0.00, ...
                'bottomMargin',   0.09, ...
                'topMargin',      0.02);

        case '1x2 tall'
            f.figureSize = [750 680];
            f.subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                'rowsNum', 1, ...
                'colsNum', 2, ...
                'heightMargin',  0.00, ...
                'widthMargin',    0.06, ...
                'leftMargin',     0.09, ...
                'rightMargin',    0.00, ...
                'bottomMargin',   0.10, ...
                'topMargin',      0.02);

        case '1x2 wide'
            f.figureSize = [750 450];
            f.subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                'rowsNum', 1, ...
                'colsNum', 2, ...
                'heightMargin',  0.00, ...
                'widthMargin',    0.06, ...
                'leftMargin',     0.09, ...
                'rightMargin',    0.01, ...
                'bottomMargin',   0.13, ...
                'topMargin',      0.04);

        case '1x2 ultra wide'
            f.figureSize = [1400 500];
            f.subplotPosVectors = NicePlot.getSubPlotPosVectors(...
                'rowsNum', 1, ...
                'colsNum', 2, ...
                'heightMargin',  0.00, ...
                'widthMargin',    0.06, ...
                'leftMargin',     0.09, ...
                'rightMargin',    0.01, ...
                'bottomMargin',   0.13, ...
                'topMargin',      0.04);

        otherwise
            error('unknown panelLayout: ''%s''.', panelLayout);
    end

    
end