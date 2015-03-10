function rgcDisplayParameter(rgcP,fieldName,varargin)
% Display function for the rgcParameters class
%
%   res = rgcDisplayParameter(rgcP,fieldName,varargin)
%
% rgcP:      the parameter object (from the class rgcParameters)
% fieldName: string describing the field that we want to display
%
% This function is not normally called directly.  You should generally use
% the function rgcP.Display(fieldName);
%
% List of display fields:
%   does not work with every field so far
%
% Examples:
%  rgcP = rgcParameters;
%  rgcP.Display('name')
%  rgcP.Display('dt');
%
%  rgcP.set('data',rand(100,100,5));
%  rgcP.Display('data'); % calls this function equivalent to 
%                         % rgcDisplayParameter(rgcP,'data');
%
%
% If  an alias is missing for a field,  add it in rgcMapParameterField But
% not in this function, which only uses the standard names of the object.

if notDefined('rgcP'), error('Parameter object needed'); end
if notDefined('fieldName'), error('Field name needed');   end

% we map field name to a standard field
correctfield = rgcMapParameterField(stringFormat(fieldName));

% correctfield is a standard way of calling the field
switch(correctfield)
    case {'name'}
        fprintf('rgcP name: %s\n',rgcP.name);
        
    case {'noise'}
        rgcP.noise

    case {'data'}
        data = rgcP.get('data');
        figure();imagesc(sum(data,3)/size(data,3));title('data averaged over time');
        mplay(255*uint8(data/max(data(:))));
        
    case {'sensor'}
        rgcP.sensor
        
    case {'oi'}
        rgcP.oi
        
    case {'dt'}
        if notDefined('varargin')
            varargin = {{'ms'}};
        end
        switch lower(varargin{1}{1})
            case{'ms'}
                fprintf('Temporal Sampling: %d ms\n',rgcP.dT);
            case{'s'}
                dt = rgcP.dT/1000;
                fprintf('Temporal Sampling: %d s\n',dt);
            case{'min'}
                dt = rgcP.dT/(1000*60);
                fprintf('Temporal Sampling: %d min\n',dt);
            case{'h'}
                dt = rgcP.dT/(1000*60*60);
                fprintf('Temporal Sampling: %d hour\n',dt);
            case{'d'}
                dt = rgcP.dT/(1000*60*60*24);
                fprintf('Temporal Sampling: %d day\n',dt);
            otherwise
                fprintf('Unknown time unit');
        end
        
    case{'t'}
        rgcP.t
        
    case {'trdur'}
        rgcP.trDur
        
    case{'tshift'}
        rgcP.tShift
        
    case{'linf'}
        rgcP.linF
        
    case {'conespacing'}
        if notDefined('varargin')
            varargin = {{'um'}};
        end
        switch lower(varargin{1}{1})
            case{'angstrom'}
                fprintf('Cone Spacing: %d angstrom \n',10000*rgcP.coneSpacing);
            case{'nm'}
                fprintf('Cone Spacing: %d nm \n',1000*rgcP.coneSpacing);
            case{'um'}
                fprintf('Cone Spacing: %d um \n',rgcP.coneSpacing);
            case{'mm'}
                fprintf('Cone Spacing: %d mm \n',rgcP.coneSpacing/1000);
            case{'cm'}
                fprintf('Cone Spacing: %d cm \n',rgcP.coneSpacing/10000);
            case{'m'}
                fprintf('Cone Spacing: %d m \n',rgcP.coneSpacing/1000000);
            case{'km'}
                fprintf('Cone Spacing: %d km \n',rgcP.coneSpacing/1000000000);
            otherwise
                fprintf('Unknown length unit');
        end
        
    case {'conespacingdefault'}
        fprintf('Default Cone Spacing: %d um \nUsed only if no cones spacing is specified in the data.',rgcP.coneSpacingDefault);
        
    case {'cellspacing'}
        if notDefined('varargin')
            varargin = {{'um'}};
        end
        nL = rgcP.nLayer;
        cs = rgcP.cellSpacing;
        for ii = 1:nL
            fprintf('Layer %d\n',ii);
            switch lower(varargin{1}{1})
                case{'angstrom'}
                    fprintf('Cell Spacing: %d angstrom \n',10000*cs{ii});
                case{'nm'}
                    fprintf('Cell Spacing: %d nm \n',1000*cs{ii});
                case{'um'}
                    fprintf('Cell Spacing: %d um \n',cs{ii});
                case{'mm'}
                    fprintf('Cell Spacing: %d mm \n',cs{ii}/1000);
                case{'cm'}
                    fprintf('Cell Spacing: %d cm \n',cs{ii}/10000);
                case{'m'}
                    fprintf('Cell Spacing: %d  m \n',cs{ii}/1000000);
                case{'km'}
                    fprintf('Cell Spacing: %d km \n',cs{ii}/1000000000);
                otherwise
                    fprintf('Unknown length unit');
            end
        end
        
    case {'actthreshparam'}
        atp = rgcP.get('actthreshparam');
        for ii = 1:rgcP.get('nlayers')
            fprintf('Layer %d. The activation threshold parameter is: %d\n',ii,atp{ii});
        end
        
    case {'meanv'}
        rgcP.meanV           

    case {'stdv'}
        rgcP.stdV            

    case {'savedir'}
        rgcP.saveDir
        
    case {'cutoff'}
        fprintf('The cutoff parameter is set to %d ',rgcP.cutoff);
        fprintf('This means that any connection below %d is cut\n',rgcP.cutoff);
    
    case {'distancefunction'}
        rgcP.distanceFunction

    case {'wrange'}
        rgcP.wRange
        
    case {'nrange'}
        rgcP.nRange

    case{'layer'}
        nL = rgcP.get('nLayers');
        for ii = 1:nL
            if isequal(lower(rgcP.getLayer(ii).type),'on')
                fprintf('%d-th layer is an ON layer \n',ii);
            else
                fprintf('%d-th layer is an OFF layer \n',ii);
            end
        end
        if (nL==0)
            fprintf('No layer \n');
        end
        
    case {'nframe'}
        fprintf('There are %d frames in the data \n',rgcP.nFrame);
        
    case {'imgidx'}
        rgcP.imgIdx

    case {'imgdur'}
        rgcP.imgDur 

    case {'cellgrid'}
        cg = rgcP.cellGrid;
        [a b] = meshgrid(cg{1},cg{2});
        scatter(a(:),b(:),'x');
        title('Cell grid expressed in um');
        
    case {'conegrid'}
        cg = rgcP.coneGrid;
        [a b] = meshgrid(cg{1},cg{2});
        scatter(a(:),b(:),'x');
        title('Cone grid expressed in um');
        
    case {'cellloc'}
        cl = rgcP.cellLoc;
        scatter(cl(:,1),cl(:,2),'x');
        title('Cell grid expressed in um');
        
    case {'gridsize'}
        fprintf('The rgc grid size is [%d %d] cells\n',rgcP.gridSize);
        
    case {'nlayers'}
        % rgcP.Display('nLayers')
        nL = rgcP.get('nLayers');
        if (nL == 0)
            fprintf('No layer \n');
        elseif (nL == 1)
            fprintf('There is one layer in the network \n');
        else
            fprintf('There are %d layers in the network \n',nL);
        end
        
    case {'rf'}
        % rgcP.Display('RF')
        % layerNumber = 1; rgcP.Display('RF',layerNumber);
        % rgcP.Display('RF',[],'sum')
        % rgcP.Display('RF',1,'surround');
        nL = rgcP.get('nLayers');
        if (nL==0)
            error('There is no layer, therefore no RF to display.');
        end
        
        % Parse the arguments
        if isempty(varargin) 
            thisLayer = 1:nL; whichComponent = 'sum';
        else
            v = varargin{1};
            if isempty(v{1}), thisLayer = 1:nL;
            else              thisLayer = v{1};
            end
            if length(v) > 1, whichComponent = v{2};
            else              whichComponent = 'sum';
            end
        end

        % Do this plotting
        cs = rgcP.get('coneSpacing');
        for ii=thisLayer
            RF = rgcP.getLayer(ii).get('RF');
            sz = [size(RF,2) size(RF,1)];
            X = 1:sz(1); X = cs*(X-mean(X));
            Y = 1:sz(2); Y = cs*(Y-mean(Y));
            switch(lower(whichComponent))
                case {'rf1','center'}
                    figure;mesh(X,Y,RF(:,:,1));
                    title('first RF component');
                    xlabel('um');
                    ylabel('um');
                case {'rf2','surround'}
                    figure;mesh(X,Y,RF(:,:,2));
                    title('second RF component');
                    xlabel('um');
                    ylabel('um');
                case 'sum'
                    figure;mesh(X,Y,sum(RF,3));
                    title('sum of the RF components');
                    xlabel('um');
                    ylabel('um');
                    figure;imagesc(X,Y,sum(RF,3));
                    title('sum of the RF components');
                    xlabel('um');
                    ylabel('um');
                otherwise
            end
        end
        
    case {'intr'}
        % rgcP.Display('inTR')
        % layerNumber = 1; rgcP.Display('RF',layerNumber);
        nL = rgcP.get('nLayers');
        if (nL==0)
            error('There is no layer, therefore no RF to display.');
        end
        
        % Parse the arguments
        if isempty(varargin) 
            thisLayer = 1:nL;
        else
            thisLayer = varargin{1}{1};
        end

        % Do this plotting
        for ii=thisLayer
            rgcP.getLayer(ii).Display('inTR')
        end
        
    
    case {'impr'}
        nL = rgcP.get('nLayers');
        for ii = 1:nL
            impR = rgcP.getLayer(ii).get('inTR');
            figure;plot(impR(:,1),'r');hold on;
            plot(impR(:,2),'b');hold off;
            legend('First Component','Second Component');
            title('Impulse Response Functions');
        end
        
    case {'cfa'}
        if isempty(rgcP.sensor)
            error('There is no sensor associated with the data');
        else
            figure;imagesc(rgcP.sensor.cfa.pattern);
        end
        
    case{'eyemovetype'}
        rgcP.eyeMoveType
        
    case {'maxeyemove'}
        rgcP.maxEyeMove
    
    case {'eyemovestd'}
        rgcP.eyeMoveStd
    
    case {'noeyemove'}
        rgcP.noEyeMove
    
    case {'eyepath'}
        rgcP.eyePath
        
    otherwise
        nL = rgcP.get('nLayers');
        for ii = 1:nL
            fprintf('layer %d:',ii);fprintf('\n');
            rgcP.getLayer(ii).Display(correctfield);
        end
        
end