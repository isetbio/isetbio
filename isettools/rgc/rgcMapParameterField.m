function res = rgcMapParameterField(fieldName)
% Maps fieldName to a standard format, implementing aliases 
% 
%    stdField = rgcMapParameterField(fieldName);
%
% Add aliases for both rgc and layer parameters in this function.
%
% The standard format is lower case with no spaces.  
%
% By using this function, we can refer to parameters in clearer text. For
% example, we can use 'Starting Noise' to indicate the parameter noise. We
% can use 'simulationData' or 'Simulation Data' to indicate the data field.
% 
%
% Examples:
%    stdField = rgcMapParameterField('Spatial Sampling')
%    rgcMapParameterField('Name')
%    rgcMapParameterField('Temporal Sampling Rate')
%    rgcMapParameterField('Connection Matrix')
%
% (c) Synapse project IBM/Darpa 2009

fieldName = stringFormat(fieldName);

switch(fieldName)
    case {'name','celltype'}
        res = 'name';
        
    case {'noise','startingnoise','noiseatthebeginning'}
        res = 'noise';

    case {'absorptions','coneabsorptions'}
        res = 'absorptions';
        
    case {'cvolts','conevolts','conevoltages'}
        res = 'cvolts';
        
    case {'sensor'}
        res = 'sensor';
        
    case {'oi','opticalimage'}
        res = 'oi';
        
    case {'dt','temporalsampling','temporalsamplingrate'}
        res = 'dt';
        
    case {'dtdefault','dt_default', 'defaultdt', 'default_dt'}
        res = 'dT_default';
        
    case {'conespacing','coneaperture'}
        res = 'conespacing';
        
    case {'conespacingdefault','defaultconespacing'}
        res = 'conespacingdefault';
        
    case {'cellspacing','rgcspacing'}
        res = 'cellspacing';
        
    case {'fps','framespersecond','framepersecond','framerate'}
        res = 'fps';
        
    case {'vswing','voltageswing'}
        % Range is 0 to voltage swing for the RGC in this layer
        res = 'vswing';
        
    case {'rgcvoltthresh','rgcvt'}
        % Spiking threshold in volts required for a spike
        res = 'rgcvoltthresh';
        
    case {'meanv','randombasecurrentmean','meanofrandombasecurrent','meanoftherandombasecurrent'}
        res = 'meanv';           

    case {'stdv','sddv','randombasecurrentstandarddeviation','randombasecurrentstd','randombasecurrentsdd'}
        res = 'stdv';            

    case {'savedir','savedirectory','savingdirectory'}
        res = 'savedir';
        
    case {'cutoff','networkcutoff'}
        res = 'cutoff';
    
    case {'distancefunction','distfunc','networkdistancefunction','netdistfunc'}
        res = 'distancefunction';
        
    case {'wrange'}
        res = 'wrange';
        
    case {'nrange'}
        res = 'nrange';

    case{'layer','layers'}
        res = 'layer';
        
    case {'nframe','nframes','numberofframes'}
        % Number of temporal frames?
        res = 'nframe';
        
    case {'gridsize'}
        res = 'gridsize';
        
    case {'cellgrid','cellmesh','cellsgrid','cellsmesh'}
        res = 'cellgrid';
        
    case {'conegrid','conemesh','conesgrid','conesmesh'}
        res = 'conegrid';
        
    case {'nlayers','nlayer','numberoflayers','nl'}
        res = 'nlayers';

    case{'type'}
        res = 'type';
        
    case{'parent'}
        res = 'parent';
        
    case{'rfstd','rfsstd','rfstandarddeviation','rfstdparam','rfstdparams'} % standard deviation sx^2 sy^2 rotationAngle
        res = 'rfstd';
        
    case{'rfcovariance','rfcovariancematrix','rfcov'} % juste the covariance matrix
        res = 'rfcov';
        
    case{'rfscoeffs','rfcoeffs','rfscoeff','rfcoeff','rfcoefficients','rfcoefs'}
        res = 'rfscoeffs';
    
    case{'cellloc','celllocation','celllocations'}
        res = 'cellloc';
        
    case{'trdur'}
        res = 'trdur';
        
    case {'trsamples','timesamples','temporalresponsetimes'}
        res = 'trsamples';
        
    case{'centertr'}
        res = 'centertr';
        
    case{'centerbdiv'}
        res = 'centerbdiv';
        
    case{'centerf'}
        res = 'centerF';
        
    case{'centernormf'}
        res = 'centernormf';
        
    case{'centergamma'}
        res = 'centergamma';
        
    case{'surroundtr'}
        res = 'surroundtr';
        
    case{'surroundbdiv'}
        res = 'surroundbdiv';
        
    case{'surroundf'}
        res = 'surroundf';
        
    case{'surroundnormf'}
        res = 'surroundnormf';
        
    case{'surroundgamma'}
        res = 'surroundgamma';
        
    case{'gammaftr'}
        res = 'gammaftr';
        
    case{'gammabdiv'}
        res = 'gammabdiv';
        
    case{'gammaf'}
        res = 'gammaf';
        
    case{'gammanormf'}
        res = 'gammanormf';
        
    case{'gammagamma'}
        res = 'gammagamma';
                       
    case{'tshift'}
        res = 'tshift';
        
    case{'linf'}
        res = 'linf';
        
    case{'intr','inputtemporalresponse','temporalinputresponse'}
        res = 'intr';
        
    case{'fbtr','feedbacktemporalresponse'}
        res = 'fbtr';
        
    case{'couplingtemporalresponse','cptr'}
        res = 'cptr';
        
    case {'numberofrgc','nrgc','numberrgc','numberofcell','nrgcs','numberrgcs','numberofrgcs','numberofcells'}
        res = 'nrgc';
    
    case {'rfcomponents','rfparts','rfcentersurround','rf'}
        % A 3D matrix of center and surround separated by components
        % rf(:,:,ii) for the ith component.
        res = 'rfcomponents';
        
    case{'rfsum','receptivefield'}
        % The sum of the rf components, rf(:,:,ii), sum(rf,3)
        res = 'rfsum';
        
    case {'irgc','indexrgc','indicesrgc','indexofcells','indicesofcells'}
        res = 'irgc';
        
    case {'impr','imprf','imprfs','impulseresponsefunction','impulseresponsefunctions'}
        res = 'impr';
        
    case {'noiseframes','noiseframe'}
        res = 'noiseframes';
        
    case {'cfa','cfagrid','cfapattern'}
        res = 'cfa';
        
    case {'layercenter','layercenterpos','layercenterposition'}
        res = 'layercenter';
       
    case {'rgcts','rgcvoltagetimeseries','rgcvolttimeseries'}
        res = 'rgcvtimeseries';
               
    case {'spiketimes','spiketimeseries','spikes','currentspiketimes','currentspkts','currentspktsstate','spktsstate'}    
        res = 'currentspkts';
        
    case {'lineartimeseries','lints','linearts','currentlints','currentlintsstate','lintsstate'}    
        res = 'currentlints';
        
    case {'componentts','componenttimeseries','rfcomponentts','rfcomponenttimesedries'}    
        res = 'rfcomponentts';

    case {'currentconnec','currentconnecstate','connecstate',...
            'currentconnect','currentconnectstate','connectstate'...
            'currentconnecmatrix','currentconnecmatrixstate','connecmatrixstate',...
            'connectionmatrix','connectionm',...
            'currentconnectmatrix','currentconnectmatrixstate','connectmatrixstate'}    
        res = 'currentconnec';
        
    case {'bordersize','sizeofborder'}
        res = 'bordersize';
        
    case {'overridesizedefault'}
        res = 'overridesizedefault';
        
    case {'overridensize'}
        res = 'overridensize';
        
    case{'eyemovetype','eyemovementtype','typeofeyemovement'}
        res = 'eyemovetype';
        
    case {'maxeyemove','maximumeyemovement'}
        res = 'maxeyemove';
    
    case {'eyemovestd','eyemovementstd','eyemovementstandartdeviation'}
        res = 'eyemovestd';
    
    case {'noeyemove','noeyemovement'}
        res = 'noeyemove';
    
    case {'eyepath','eyemovepath','eyemovementpath'}
        res = 'eyepath';
        
    otherwise
        res = fieldName;
end

% Make sure there are no caps in the return, and no spaces
res = stringFormat(res);

return
