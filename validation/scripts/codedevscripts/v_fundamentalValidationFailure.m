function varargout = v_fundamentalValidationFailure(varargin)
%
% Example validation script that demonstrates usage of the fundemantal failure feature. 
%

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)
      
    % Simulate fundamental failure here
    if (true)
        UnitTest.validationRecord('FUNDAMENTAL_CHECK_FAILED', 'Do not panic. This test is working if you see a big red fundamental failure message here.');
        return;
    end
    
    UnitTest.validationRecord('PASSED',  'all right to here');
    UnitTest.validationData('dummyData', ones(100,10));
    
    % Plotting
    if (runTimeParams.generatePlots)
       figure(1);
       clf;
       plot(1:10, 1:10, 'r-');
       axis 'square'
       drawnow;
    end
    
end