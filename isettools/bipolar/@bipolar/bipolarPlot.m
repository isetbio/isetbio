function hdl = bipolarPlot(obj, pType,varargin)
% Plot the values from the bipolar object
% 
%    hdl = bp.plot(parameter)
%
% Plot types
%   response center
%   response surround
%   response
%   movie response 
%   Time series
% 
% Examples:
%
% 
% 5/2016 JRG (c) isetbio team

%% Parse inputs

p = inputParser; 
p.CaseSensitive = false; 
p.FunctionName = mfilename;
p.KeepUnmatched = true;

% Make key properties that can be set required arguments, and require
% values along with key names.
allowableFields = {...
    'response',...
    'responseCenter','bipolarresponsecenter',...    
    'responseSurround','bipolarresponsesurround'...
    'movieresponse'
    };
p.addRequired('pType',@(x) any(validatestring(ieParamFormat(x),allowableFields)));

% Parse pType
p.parse(pType,varargin{:}); 

%% Create window
hdl = vcNewGraphWin([],'upperLeftBig');

sz = size(obj.responseCenter);

% Programming:
% We need to get the units of time from the object, not as per below.

% Options
switch ieParamFormat(pType);  
    case{'responsecenter'}
        
        responseRS = reshape(obj.responseCenter,sz(1)*sz(2),sz(3));
        plot(.001*(1:sz(3)),responseRS);
        xlabel('Time (sec)');
        ylabel('Response (AU)');
        title('Bipolar Mosaic Response');
        
    case{'responsesurround'}
        
        responseRS = reshape(obj.responseSurround,sz(1)*sz(2),sz(3));
        plot(.001*(1:sz(3)),responseRS);
        xlabel('Time (sec)');
        ylabel('Response (AU)');
        title('Bipolar Mosaic Response');
        
    case{'response'}
        
        response = reshape(obj.get('response'),sz(1)*sz(2),sz(3));
        plot(.001*(1:sz(3)),response);
        xlabel('Time (sec)');
        ylabel('Response (AU)');
        title('Bipolar Mosaic Response');
        
    case{'movieresponse'}
        % Pass the varargin along
        if ~isempty(varargin) && length(varargin) == 1
            % Params are coded in a single struct
            varargin{1}.hf = hdl;
            ieMovie(obj.get('response'),varargin{:});
        else
            % List of params
            ieMovie(obj.get('response'),'hf',hf,varargin{:});
        end
end

% set(gca,'fontsize',16);

end
