function [osCurrent, obj] = coneAdaptAlt(sensor, typeAdapt)

if strcmpi('linear',typeAdapt)
    
    linearOS = osLinear('noiseFlag', 1); %osLinear
       
    % Compute linear outer segment response.
    linearOS = osLinearCompute(linearOS, sensor);
    % params.offset = 0;
    % linearOS = osLinearCompute(linearOS, sensor, params);
    osLinearGet(linearOS, 'noiseFlag');
    
    % Plot results.
    osLinearPlot(linearOS, sensor);
    
    obj = linearOS;
    osCurrent = obj.ConeCurrentSignal;
    
elseif strcmpi('rieke', typeAdapt) || strcmpi('biophys', typeAdapt)
    
    adaptedOS = osBioPhys(); % osBioPhys
    adaptedOS = osBioPhysSet(adaptedOS, 'noiseFlag', 1);
    
    % Compute nonlinear outer segment response.
    adaptedOS = osBioPhysCompute(adaptedOS, sensor);
    % params.bgVolts = 0; params.offset = 0;
    % adaptedOS = osBioPhysCompute(adaptedOS, sensor, params);
    osBioPhysGet(adaptedOS, 'noiseFlag');
    
    % Plot results
    osBioPhysPlot(adaptedOS, sensor);
    
    obj = adaptedOS;
    osCurrent = obj.ConeCurrentSignal;
    
else
    
    fprintf('error: typeAdapt string must be either ''linear'' or ''rieke'' ');
    
end
