% AGGREGATE analyzed data from different runs at multiple eccentricities
function [surroundOptimizationStrategyToBeAggregated, ...
     targetNumberOfMappedCellsToBeAggregated, ...
     targetHorizontalEccentricitiesToBeAggregated, ...
     targetMosaicSizesToBeAggregated, ...
     targetMappedPositionDegsToBeAggregated ] = assembleEccentricityMatrices()
    

    surroundOptimizationStrategyToBeAggregated = { ...
        'LowH1paramsNarrowVisualSTFparamTolerance', ...  %  @  0.0
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -0.1
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -0.2
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -0.3
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -0.4
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -0.5
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -0.6
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -0.8
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -0.95
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -1.1
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -1.25
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -1.4
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -1.55
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -1.7
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -1.85
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -2.0
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -2.25
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -2.5
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -2.75
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -3.0
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -3.25
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -3.5
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -3.75
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -4.0
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -4.5
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -5.0
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -5.5
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -6.0
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -6.5
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -7.0
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -7.5
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -8.0
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -8.5
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -9.0
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -9.5
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -10.0
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -11.0
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -11.4
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -11.75
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -12.5
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -13.5
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -14.0
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -15.0
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -15.6
        'LowH1paramsNarrowVisualSTFparamTolerance' ...   %  @ -16.25
        'MidH1paramsNarrowVisualSTFparamTolerance', ...  %  @ -16.25
        'MidH1paramsNarrowVisualSTFparamTolerance', ...  %  @ -17.5
        'MidH1paramsNarrowVisualSTFparamTolerance', ...  %  @ -18.75
        'MidH1paramsNarrowVisualSTFparamTolerance', ...  %  @ -20.0
        'MidH1paramsNarrowVisualSTFparamTolerance', ...  %  @ -21.25
        'UpperH1paramsNarrowVisualSTFparamTolerance' ... %  @ -22.5 
        'UpperH1paramsNarrowVisualSTFparamTolerance' ... %  @ -23.75 
        'UpperH1paramsNarrowVisualSTFparamTolerance' ... %  @ -25.0 
        'UpperH1paramsNarrowVisualSTFparamTolerance' ... %  @ -26.25 
        'UpperH1paramsNarrowVisualSTFparamTolerance' ... %  @ -27.5 
        'UpperH1paramsNarrowVisualSTFparamTolerance' ... %  @ -29.25
        'UpperH1paramsNarrowVisualSTFparamTolerance' ... %  @ -32.0 
        };

    targetNumberOfMappedCellsToBeAggregated = { ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        144 ...
        };

    targetHorizontalEccentricitiesToBeAggregated = [...
          0 ... %  @   @ 0.0
          0 ... %  @  -0.1
          0 ... %  @  -0.2
          0 ... %  @  -0.3
          0 ... %  @  -0.4
          0 ... %  @  -0.5
          0 ... %  @  -0.6
          0 ... %  @  -0.8
          0 ... %  @  -0.95
         -2.0 ... %  @  -1.1
         -2.0 ... %  @  -1.25
         -2.0 ....%  @  -1.4
         -2.0 ... % -1.55
         -2.0 ....%  @  -1.7
         -2.0 ....%  @  -1.85
         -2.0 ....%  @  -2.0
         -2.0 ... %  @  -2.25
         -2.0 ... %  @  -2.5
         -2.0 ... %  @  -2.75
         -4.0 ... %  @  -3.0
         -4.0 ... %  @  -3.25
         -4.0 ... %  @  -3.5
         -4.0 ... %  @  -3.75
         -4.0 ... %  @  -4.0
         -4.0 ... %  @  -4.5
         -4.0 ... %  @  -4.0
         -7.0 ... %  @  -5.5
         -7.0 ... %  @  -6.0
         -7.0 ... %  @  -6.5
         -7.0 ... %  @  -7.0
         -7.0 ... %  @  -7.5
         -7.0 ... %  @  -8.0
         -7.0 ... %  @  -8.5
        -10.0 ... %  @  -9.0
        -10.0 ... %  @  -9.5
        -10.0 ... %  @  -10.0
        -10.0 ... %  @  -11.0
        -10.0 ... %  @ -11.4
        -14.0 ... %  @ -11.75
        -14.0 ... %  @ -12.5
        -14.0 ... %  @ -13.5
        -14.0 ... %  @ -14.0
        -14.0 ... %  @ -15.0
        -14.0 ... %  @ -15.6
        -14.0 ... %  @ -16.25
        -19.0 ....%  @ -16.25
        -19.0 ... %  @ -17.5
        -19.0 ....%  @ -18.75
        -19.0 ... %  @ -20.0
        -19.0 ... %  @ -21.25
        -25.0 ... %  @ -22.5
        -25.0 ... %  @ -23.75
        -25.0 ... %  @ -25.0
        -25.0 ... %  @ -26.25
        -25.0 ... %  @ -27.5
        -32.0 ... %  @ -29.25
        -32.0 ... %  @ -32.0
    ];

    targetMosaicSizesToBeAggregated = 1 + [...
        1 1; ...  %  @  0.0
        1 1; ...  %  @ -0.1
        1 1; ...  %  @ -0.2
        1 1; ...  %  @ -0.3
        1 1; ...  %  @ -0.4
        1 1; ...  %  @ -0.5
        1 1; ...  %  @ -0.6
        1 1; ...  %  @ -0.8
        1 1; ...  %  @ -0.95
        1 1; ...  %  @ -1.1
        1 1; ...  %  @ -1.25
        1 1; ...  %  @ -1.4
        1 1; ...  % -1.55
        1 1; ...  %  @ -1.7
        1 1; ...  %  @ -1.85
        1 1; ...  %  @ -2.0
        1 1; ...  %  @ -2.25
        1 1; ...  %  @ -2.5
        1 1; ...  %  @ -2.75
        2 2; ...  %  @ -3.0
        2 2; ...  %  @ -3.25
        2 2; ...  %  @ -3.5
        2 2; ...  %  @ -3.75
        2 2; ...  %  @ -4.0
        2 2; ...  %  @ -4.5
        2 2; ...  %  @ -5.0
        3 3; ...  %  @ -5.5
        3 3; ...  %  @ -6.0
        3 3; ...  %  @ -6.5
        3 3; ...  %  @ -7.0
        3 3; ...  %  @ -7.5
        3 3; ...  %  @ -8.0
        3 3; ...  %  @ -8.5
        3 3; ...  % @  -9.0
        3 3; ...  % @  -9.5
        3 3; ...  % @ -10.0
        3 3; ...  % @ -11.0
        3 3; ...  % @ 11.4
        4 4; ...  % @ -11.75
        4 4; ...  % @ -12.5
        4 4; ...  % @ -13.5
        4 4; ...  % @ -14.0
        4 4; ...  % @ -15.0
        4 4; ...  % @ -15.6
        4 4; .... % @ -16.25
        5 5; .... % @ -16.25
        5 5; .... % @ -17.5
        5 5; .... % @ -18.75
        5 5; .... % @ -20.0
        5 5; ...  % @ -21.25
        6 6; ...  % @ -22.5
        6 6; ...  % @ -23.75
        6 6; ...  % @ -25.0
        6 6; ...  % @ -26.75
        6 6; ...  % @ -27.5
        8 8; ...  % @ -29.25
        8 8 ...   % @ -32.0
    ];

    % Mapped positions (by default the horizontal centers of the mosaics, with vertical ecc = 0)
    targetMappedPositionDegsToBeAggregated(:,1) = targetHorizontalEccentricitiesToBeAggregated;
    targetMappedPositionDegsToBeAggregated(:,2) = zeros(numel(targetHorizontalEccentricitiesToBeAggregated),1);

    % Modified mapped positions, other than the default positions, which are the mosaic centers
    targetMappedPositionDegsToBeAggregated(1,1)  =  0.0;     % within the 0 deg mosaic 
    targetMappedPositionDegsToBeAggregated(2,1)  = -0.1;     % within the 0 deg mosaic
    targetMappedPositionDegsToBeAggregated(3,1)  = -0.2;     % within the 0 deg mosaic
    targetMappedPositionDegsToBeAggregated(4,1)  = -0.3;     % within the 0 deg mosaic
    targetMappedPositionDegsToBeAggregated(5,1)  = -0.4;     % within the 0 deg mosaic
    targetMappedPositionDegsToBeAggregated(6,1)  = -0.5;     % within the 0 deg mosaic
    targetMappedPositionDegsToBeAggregated(7,1)  = -0.6;     % within the 0 deg mosaic
    targetMappedPositionDegsToBeAggregated(8,1)  = -0.8;     % within the 0 deg mosaic
    targetMappedPositionDegsToBeAggregated(9,1)  = -0.95;    % within the 0 deg mosaic
    
    targetMappedPositionDegsToBeAggregated(10,1)  = -1.1;    % within the 2 deg mosaic
    targetMappedPositionDegsToBeAggregated(11,1)  = -1.25;   % within the 2 deg mosaic 
    targetMappedPositionDegsToBeAggregated(12,1)  = -1.4;    % within the 2 deg mosaic
    targetMappedPositionDegsToBeAggregated(13,1)  = -1.55;   % within the 2 deg mosaic
    targetMappedPositionDegsToBeAggregated(14,1) =  -1.7;    % within the 2 deg mosaic
    targetMappedPositionDegsToBeAggregated(15,1) =  -1.85;   % within the 2 deg mosaic
    targetMappedPositionDegsToBeAggregated(16,1) =  -2.0;    % within the 2 deg mosaic
    targetMappedPositionDegsToBeAggregated(17,1) =  -2.25;   % within the 2 deg mosaic
    targetMappedPositionDegsToBeAggregated(18,1) =  -2.5;    % within the 2 deg mosaic
    targetMappedPositionDegsToBeAggregated(19,1) =  -2.75;   % within the 2 deg mosaic
    
    targetMappedPositionDegsToBeAggregated(20,1) =  -3.0;    % within the 4 deg mosaic 
    targetMappedPositionDegsToBeAggregated(21,1) =  -3.25;   % within the 4 deg mosaic 
    targetMappedPositionDegsToBeAggregated(22,1) =  -3.5;    % within the 4 deg mosaic 
    targetMappedPositionDegsToBeAggregated(23,1) =  -3.75;   % within the 4 deg mosaic 
    targetMappedPositionDegsToBeAggregated(24,1) =  -4.0;    % within the 4 deg mosaic
    targetMappedPositionDegsToBeAggregated(25,1) =  -4.5;    % within the 4 deg mosaic 
    targetMappedPositionDegsToBeAggregated(26,1) =  -5.0;    % within the 4 deg mosaic
    
    targetMappedPositionDegsToBeAggregated(27,1) =  -5.5;    % within the 7 deg mosaic
    targetMappedPositionDegsToBeAggregated(28,1) =  -6.0;    % within the 7 deg mosaic
    targetMappedPositionDegsToBeAggregated(29,1) =  -6.5;    % within the 7 deg mosaic
    targetMappedPositionDegsToBeAggregated(30,1) =  -7.0;    % within the 7 deg mosaic
    targetMappedPositionDegsToBeAggregated(31,1) =  -7.5;    % within the 7 deg mosaic 
    targetMappedPositionDegsToBeAggregated(32,1) =  -8.0;    % within the 7 deg mosaic 
    targetMappedPositionDegsToBeAggregated(33,1) =  -8.5;    % within the 7 deg mosaic

    targetMappedPositionDegsToBeAggregated(34,1) =  -9.0;    % within the 10 deg mosaic
    targetMappedPositionDegsToBeAggregated(35,1) =  -9.5;    % within the 10 deg mosaic 2.5mm
    targetMappedPositionDegsToBeAggregated(36,1) = -10.0;    % within the 10 deg mosaic 2.7 mm
    targetMappedPositionDegsToBeAggregated(37,1) = -11.0;    % within the 10 deg mosaic  -3.0 mm - not very uniform rfs
    targetMappedPositionDegsToBeAggregated(38,1) = -11.4;    % ithin the 10 deg mosaic  -3.15 mm

    targetMappedPositionDegsToBeAggregated(39,1) = -11.75;   % within the -14 deg mosaic  -3.3mm
    targetMappedPositionDegsToBeAggregated(40,1) = -12.5;    % within the -14 deg mosaic  -3.4mm
    targetMappedPositionDegsToBeAggregated(41,1) = -13.5;    % within the -14 deg mosaic  -3.7mm
    targetMappedPositionDegsToBeAggregated(42,1) = -14.0;    % within the -14 deg mosaic  -3.8 mm   
    targetMappedPositionDegsToBeAggregated(43,1) = -15.0;    % within the -14 deg mosaic  -4.1 mm
    targetMappedPositionDegsToBeAggregated(44,1) = -15.6;    % within the -14 deg mosaic  -4.3 mm
    targetMappedPositionDegsToBeAggregated(45,1) = -16.25;   % within the -14 deg mosaic  -4.4 mm

    targetMappedPositionDegsToBeAggregated(46,1) = -16.25;   % within the -19 deg mosaic  -4.6 mm
    targetMappedPositionDegsToBeAggregated(47,1) = -17.5;    % within the -19 deg mosaic  -4.8 mm
    targetMappedPositionDegsToBeAggregated(48,1) = -18.75;   % within the -19 deg mosaic  -5.1 mm
    targetMappedPositionDegsToBeAggregated(49,1) = -20.0;    % within the -19 deg mosaic  -5.5 mm
    targetMappedPositionDegsToBeAggregated(50,1) = -21.25;   % within the -19 deg mosaic  -5.8 mm

    targetMappedPositionDegsToBeAggregated(51,1) = -22.5;    % within the -25 deg mosaic - 6.2 mm
    targetMappedPositionDegsToBeAggregated(52,1) = -23.75;   % within the -25 deg mosaic - 6.5 mm
    targetMappedPositionDegsToBeAggregated(53,1) = -25.0;    % within the -25 deg mosaic - 6.8 mm
    targetMappedPositionDegsToBeAggregated(54,1) = -26.25;   % within the -25 deg mosaic - 7.2 mm
    targetMappedPositionDegsToBeAggregated(55,1) = -27.5;    % within the -25 deg mosaic - 7.5 mm

    targetMappedPositionDegsToBeAggregated(56,1) = -29.25;   % within the -32 deg mosaic - 8.0 mm
    targetMappedPositionDegsToBeAggregated(57,1) = -32.0;    % within the -32 deg mosaic


    assert(size(targetMappedPositionDegsToBeAggregated,1) == size(targetMosaicSizesToBeAggregated,1), ...
        '''targetMappedPositionDegsToBeAggregated'' and ''targetMosaicSizesToBeAggregated'' do not match in size.');

    assert(size(targetMappedPositionDegsToBeAggregated,1) == numel(targetHorizontalEccentricitiesToBeAggregated), ...
        '''targetMappedPositionDegsToBeAggregated'' and ''targetHorizontalEccentricitiesToBeAggregated'' do not match in size.');

    assert(size(targetMappedPositionDegsToBeAggregated,1) == numel(surroundOptimizationStrategyToBeAggregated), ...
        '''targetMappedPositionDegsToBeAggregated'' and ''surroundOptimizationStrategyToBeAggregated'' do not match in size.');

    assert(size(targetMappedPositionDegsToBeAggregated,1) == numel(targetNumberOfMappedCellsToBeAggregated), ...
        '''targetMappedPositionDegsToBeAggregated'' and ''targetNumberOfMappedCellsToBeAggregated'' do not match in size.');


    if (1==2)
        % Remove some conditions
        idx = 1:numel(targetHorizontalEccentricitiesToBeAggregated);
        idx = setdiff(idx, [43]);

        surroundOptimizationStrategyToBeAggregated = surroundOptimizationStrategyToBeAggregated(idx);
        targetHorizontalEccentricitiesToBeAggregated = targetHorizontalEccentricitiesToBeAggregated(idx);
        targetNumberOfMappedCellsToBeAggregated =  targetNumberOfMappedCellsToBeAggregated(idx);
        targetMappedPositionDegsToBeAggregated = targetMappedPositionDegsToBeAggregated(idx,:);
        targetMosaicSizesToBeAggregated = targetMosaicSizesToBeAggregated(idx,:);
    end

end
