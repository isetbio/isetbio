function layerDisplay(layer,fieldName)
% implements the Display function for the rgcLayer class
%
% Obsolete?
%
% layerDisplay(layer,fieldName,value)
%
% layer: the parameter object (from the class rgcParameters)
% fieldName: string describing the field that we want
% value: value to be set to
%
% This function is equivalent to
% layer.Display(fieldName);
% which is probably the way it should be called
%
% Ex:
%
% If you feel like an alias is missing for a field, feel free to add it in rgcMapParameterField if
% it does not already exist for another field. But not in this function,
% this only uses the standard name.
%
% This function needs to be developped to display the receptive fields,
% parameters...

if notDefined('layer')
    error('layer object needed');
end

if notDefined('fieldName')
    error('Field name needed');
end

% we map field name to a standard field
correctfield = rgcMapParameterField(stringFormat(fieldName));

% correctfield is a standard way of calling the field
switch(correctfield)
    case{'name'}
        fprintf(strcat(strcat('The name of the layer is: ',layer.name),'\n'));
        
    case{'parent'}
        fprintf('The parent rgcParameters object for this layer is :\n');
        layer.parent
        
    case{'cellspacing'}
        fprintf('The cell spacing for this layer is %.3f um\n',layer.cellSpacing);
        %layer.cellSpacing
        
    case {'gridsize'}
        fprintf('The gridsize is %d by %d',layer.gridSize(1), layer.gridSize(2));
        %layer.gridSize
        
    case {'cellgrid'}
        xaxis = layer.cellGrid{1}
        yaxis = layer.cellGrid{2}
        fprintf('The grid of the locations of the cells is : \n');
        fprintf('x axis : from %.2f to %.2f um, with a step of %.2f um\n', min(xaxis), max(xaxis), xaxis(2)-xaxis(1));
        fprintf('y axis : from %.2f to %.2f um, with a step of %.2f um\n', min(yaxis), max(yaxis), yaxis(2)-yaxis(1));
        display = input('Do you want to display the axes? (1 for yes, 0 for no)');
        if display==1
           fprintf('x axis :\n');
           xaxis
           fprintf('y axis :\n');
           yaxis
        else if display ~=0
                error('input has to be 0 or 1');
            end
        end
        
    case{'cellloc'}
        % display the first ten, then ask for the rest
        fprintf('The first 10 cell locations (x,y) are : \n');
        layer.cellLoc(1:10,:)
        fprintf('There are %d cells. ',size(layer.cellLoc,1));
        display = input('Display all? (1 for yes, 0 for no)');
        if display==1
            layer.cellLoc(11:end,:)
        else if display ~=0
                error('input has to be 0 or 1');
            end
        end
        
    case{'rfcov'}
        fprintf('the covariance matrix (in um^2) is:\n');
        fprintf('for the center:\n');
        disp(layer.RFcov{1});
        fprintf('for the surround:\n');
        disp(layer.RFcov{2});
%  
%     case{'rfstd'}
%         cov = layer.RFcov;
%         lc = length(cov); % = nb components, should be 2
%         for ii = 1:lc
%             res{ii} = rgcExtractStdFromCovariance(cov{ii});
%         end
%         % assuming lc = 2
%         fprintf('for the center:\n   Standard deviations : %.2f and %.2f um \n   in the directions [%.2f %.2f] and [%.2f %.2f]\n   corresponding to an angle of : %.2f rad = %.2f deg', ...
%             res{1}.std(1),res{1}.std(2),res{1}.directions(1,1),res{1}.directions(1,2),res{1}.directions(2,1),res{1}.directions(2,2),res{1}.angle,res{1}.angle/pi*180);
%         fprintf('\n');
%         fprintf('for the surround:\n   Standard deviations : %.2f and %.2f um \n   in the directions [%.2f %.2f] and [%.2f %.2f]\n   corresponding to an angle of : %.2f rad = %.2f deg', ...
%             res{2}.std(1),res{2}.std(2),res{2}.directions(1,1),res{2}.directions(1,2),res{2}.directions(2,1),res{2}.directions(2,2),res{2}.angle,res{2}.angle/pi*180);
%         fprintf('\n');
               
    case{'rfscoeffs'}
        fprintf('The RF coefficients are %.2f for the center and %.2f for the surround.\n',layer.RFsCoeffs(1), layer.RFsCoeffs(2));
                 
    case{'trdur'}
        fprintf('Duration of the input response function, in ms:\n');
        layer.trDur
        
    case{'centertr'}
        layer.centertr
        
    case{'centerbdiv'}
        layer.centerBDiv
        
    case{'centerf'}
        layer.centerF
        
    case{'centernormf'}
        layer.centerNormF
        
    case{'centergamma'}
        layer.centerGamma
        
    case{'surroundtr'}
        layer.surroundTR
        
    case{'surroundbdiv'}
        layer.surroundBDiv
        
    case{'surroundf'}
        layer.surroundF
        
    case{'surroundnormf'}
        layer.surroundNormF
        
    case{'surroundgamma'}
        layer.surroundGamma
        
    case{'gammaftr'}
        layer.gammafTR
        
    case{'gammabdiv'}
        layer.gammaBDiv
        
    case{'gammaf'}
        layer.gammaF
        
    case{'gammanormf'}
        layer.gammaNormF
        
    case{'gammagamma'}
        layer.gammaGamma
        
    case{'dt'}
        fprintf('The layer spiking time step (dT) is %.2f ms\n',layer.dT);
        
    case{'tshift'}
        layer.tShift
        
    case{'linf'}
        layer.linF
        
    case{'intr'}
        figure;
        plot(layer.inTR);
        title('Input Temporal response');
        xlabel('ms');
        
    case {'actthreshparam'}
        layer.ActThreshParam
    
    case{'fbtr'}
        figure;
        plot(layer.fbTR);
        title('Feedback Temporal response');
        xlabel('ms');
        
    case{'cptr'}
        figure;
        plot(layer.cpTR);
        title('Coupling Temporal response');
        xlabel('ms');
        
    case{'rf'}
        RF = layer.RF;
        cs = layer.cellSpacing;
        sz = [size(RF,2) size(RF,1)];
        X = 1:sz(1); X = cs*(X-mean(X));
        Y = 1:sz(2); Y = cs*(Y-mean(Y));
        figure;mesh(X,Y,sum(RF,3));
        title('sum of the RF components');
        xlabel('um');
        ylabel('um');
        figure;imagesc(X,Y,sum(RF,3));
        title('sum of the RF components');
        xlabel('um');
        ylabel('um');
        
    case {'cutoff'}
        layer.cutoff
        
    case {'wrange'}
        layer.wRange
        
    case {'nrange'}
        layer.nRange
    
    case {'rgcvtimeseries'}
        % current voltage of each cell
        layer.rgcvTimeSeries
        
    case {'currentspikeintegral'}
        % current spike integral of the cell
        layer.currentSpikeIntegral;
        
    case {'currentspkts'}
        % spike time series, for each cell
        layer.currentSpkTS
        
    case {'currentlints'}
        % linear time series, for each cell
        layer.currentLinTS
                
    case {'currentconnec'}
        % connexion matrix
        layer.currentConnec
        
    otherwise
        error('Unknown field name');
end