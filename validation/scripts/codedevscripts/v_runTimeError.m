function varargout = v_runTimeError(varargin)
%
%  Example validation script that simulates runtime exemption. 
%

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)
    
    error('Simulating runtime error');
    
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