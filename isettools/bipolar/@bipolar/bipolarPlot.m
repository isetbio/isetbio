function bipolarPlot(obj, varargin)
% Plot the input and output of the bipolar object.
% 
% The bipolar object allows the simulated cone responses to be passed on to
% the inner retina object and to approxiately maintain its impulse
% response. This will allow us to run the nonlinear biophysical cone outer
% segment model and pass its results on to the bipolar stage and then RGCs.
% 
% 5/2016 JRG (c) isetbio team

%%
narginchk(0, Inf);
p = inputParser; p.CaseSensitive = false; p.FunctionName = mfilename;
p.KeepUnmatched = true;
% Make key properties that can be set required arguments, and require
% values along with key names.
allowableFieldsToSet = {...
    'response',...
    'responseCenter','bipolarresponsecenter',...    
    'responseSurround','bipolarresponsesurround'...
    };
p.addRequired('what',@(x) any(validatestring(ieParamFormat(x),allowableFieldsToSet)));

% Parse and put results into structure p.
p.parse(varargin{:}); params = p.Results;

vcNewGraphWin([],'upperLeftBig');

responseSize = size(obj.responseCenter);

% 
switch ieParamFormat(params.what);  
    case{'response'}
        
        responseRS = reshape(obj.responseCenter,responseSize(1)*responseSize(2),responseSize(3));

        plot(.001*(1:responseSize(3)),responseRS);
        xlabel('Time (sec)');
        ylabel('Response (AU)');
        title('Bipolar Mosaic Response');
        
    case{'responseSurround'}
        responseRS = reshape(obj.responseSurround,responseSize(1)*responseSize(2),responseSize(3));
        plot(.001*(1:responseSize(3)),responseRS);
        xlabel('Time (sec)');
        ylabel('Response (AU)');
        title('Bipolar Mosaic Response');
    case{'responseCenter'}
        responseRS = reshape(obj.responseCenter-obj.responseSurround,responseSize(1)*responseSize(2),responseSize(3));
        plot(.001*(1:responseSize(3)),responseRS);
        xlabel('Time (sec)');
        ylabel('Response (AU)');
        title('Bipolar Mosaic Response');
end
set(gca,'fontsize',16);
