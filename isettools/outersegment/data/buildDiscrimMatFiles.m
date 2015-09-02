


clear


for color_val = 1%[1 3 4]
%     for contrast = [0 0.001 0.002 0.004 0.006 0.008 0.01]% 0.05 0.1 0.15 0.2 0.4 0.6 0.8]
for contrast = [0.0002 0.0004 0.0006 0.0008]
        [color_val contrast]
        % parameters found in Fig. 6 caption
        
        params = paramsGaborColorOpponent()
        params.color_val = color_val;         % 1 = s_iso, 2 = L-M, 3 = LMS, 4 = L-M
        params.disp_movie = 0;        % display movie flag
        params.contrast = contrast;
        
        % build scene
        [scene, display] = sceneHorwitzHass(params);
        displayClose;
        
        % build optical image
        % oi  = oiCreate('wvf human');
        
        %  According to the paper, pupil size should be of area 12.6 mm^2 (4 mm
        %  pupil diameter)
        pupil_size = 4; % 4 mm diameter
        oi = oiCreate('wvf human', pupil_size);

        % build sensor
        % Gabor color-opponent present
        sensor = sensorHorwitzHass(params, scene, oi, display);   
        
        savestr = ['sensor_new_CV_' num2str(color_val) '_Cont_' sprintf('%0.4d',10000*contrast) '.mat'];
        save(savestr,'sensor','scene','oi','display','params');        
        
        % build outersegment
        noiseflag = 1;
        linearOS = osLinear('noiseFlag',noiseflag);
        % adaptedOS = osBioPhys('noiseFlag',noiseflag);
        % compute outersegment current output
        linearOS = osLinearCompute(linearOS,sensor,params);
        
        pooledData = pooledConeResponse_orig(linearOS, sensor);
        
       
        savestr2 = ['pooled_Response_new_CV_' num2str(color_val) '_Cont_' sprintf('%0.4d',10000*contrast) '.mat'];
        save(savestr2, 'pooledData');
        
        clear sensor scene oi disp params pooledData linearOS
        
        
    end
end