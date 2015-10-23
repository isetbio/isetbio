function [raster psth] = computePSTH(obj)


dt = .01; % make this a get from sensor
nCells = size(obj.cellLocation);

spikeCheck = (cellfun(@isempty,((obj.spikeResponse))));

if sum(spikeCheck(:)) ~= (length(obj.spikeResponse)*nCells(1)*nCells(2))

maxTrials = size(obj.spikeResponse,3);

for xcell = 1:nCells(1)
    for ycell = 1:nCells(2)
        
        
        for trial = 1:maxTrials
            tsp{trial} = obj.spikeResponse{xcell,ycell,trial,1};
        end
        
        
        if sum(cellfun(@isempty,tsp))~=maxTrials
            mtsp = plotraster(tsp);
        else
            mtsp = [];
        end
        
        raster{xcell,ycell} = mtsp;
        
        [psth{xcell,ycell},tt,pstv,spr] = compPSTH(mtsp*dt, .001, .002, [0 1], .005);
        
    end
end

else 
    raster = []; psth = [];

end