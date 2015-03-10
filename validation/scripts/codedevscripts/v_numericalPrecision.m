function varargout = v_numericalPrecision(varargin)
%
% Script assessing the effects of rounding at different numerical precisions.
%

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)

    
    %% Initialize ISETBIO
    s_initISET;
        
    %% Internal validations
    % The result of a computation: a vector of doubles
    gain = 10E18;
    largeResult = gain*eps*rand(1,100*1000);
    largeResult = largeResult - mean(largeResult);
    
    gain = 200;
    smallResult = gain*eps*rand(1,100*1000);
    smallResult = smallResult - mean(smallResult);
    
    decimalDigits = [10:16];
    tolerance = 500*eps;
     
    %% Conversion to N- digit precision
    for k = 1:numel(decimalDigits)
        nDigits = decimalDigits(k);
        largeResultRepresentations(k,:) = round(largeResult,nDigits);
        diff = abs(largeResult-squeeze(largeResultRepresentations(k,:)));
        if (any(diff >= tolerance ))
            UnitTest.validationRecord('SIMPLE_MESSAGE', sprintf('large result is not represented accurately with %d decimal digits.', nDigits));
        else
            UnitTest.validationRecord('SIMPLE_MESSAGE', sprintf('large result is represented accurately with %d decimal digits.', nDigits));
        end
    end
    
    for k = 1:numel(decimalDigits)
        nDigits = decimalDigits(k);
        smallResultRepresentations(k,:) = round(smallResult,nDigits);
        diff = abs(smallResult-squeeze(smallResultRepresentations(k,:)));
        if (any(diff >= tolerance ))
            UnitTest.validationRecord('SIMPLE_MESSAGE', sprintf('small result is not represented accurately with %d decimal digits.', nDigits));
        else
            UnitTest.validationRecord('SIMPLE_MESSAGE', sprintf('small result is represented accurately with %d decimal digits.', nDigits));
        end

    end
    
    %% Append to extraData
    UnitTest.extraData('largeResult', largeResult);
    UnitTest.extraData('resultRepresentations', largeResultRepresentations);
    UnitTest.extraData('smallResult', smallResult);
    UnitTest.extraData('smallRepresentations', smallResultRepresentations);
     
    %% Plotting
    if (runTimeParams.generatePlots)
        plotResults(1, decimalDigits, largeResult, largeResultRepresentations, 'Large Values');
        plotResults(2, decimalDigits, smallResult, smallResultRepresentations, 'Small Values');  
    end
    
end

%% Helper plotting function
function plotResults(figNum, decimalDigits, result, resultRepresentations, figName)
    h = figure(figNum);
    set(h, 'Position', [100 100  1310 900], 'Name', figName);
    clf;
    subplotWidth = 0.62/7;
    margin = 0.05;
    for k = 1:numel(decimalDigits)
        nDigits = decimalDigits(k);
        subplot('Position', [0.06+(k-1)*(subplotWidth+margin), 0.05, subplotWidth 0.90]);
        diff = abs(result-squeeze(resultRepresentations(k,:)));
        plot(result, diff, 'b.');
        hold on; 
        plot([min(result) max(result)], [eps eps], 'r-');
        
        hold off
        set(gca, 'XLim', [min(result) max(result)], 'YLim', [0 10000*eps], 'YScale', 'log', ...
             'YTick', eps*[1 10 100 300 1000 3000 10000], ...
             'YTickLabel', {'EPS', '10EPS', '100EPS', '300EPS', '1kEPS',  '3kEPS', '10kEPS'});
        set(gca, 'FontSize', 12, 'FontName', 'Helvetica');

        xlabel('value', 'FontSize', 14, 'FontName', 'Helvetica', 'FontWeight', 'b');
        if (mod(k-1,7) == 0)
            ylabel('| value - round(value,nDigits) |', 'FontSize', 14, 'FontName', 'Helvetica', 'FontWeight', 'b');
        end
        %   set(gca, 'YTick', []) 
        
        title(sprintf('nDigits = %d', nDigits), 'FontSize', 14, 'FontName', 'Helvetica');
    end
end
        