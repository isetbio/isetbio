% Method to compute the bandpass index of an STF as defined in 
% Lee, Shapley, Hayken, and Sun. (2012) "Spatial distributions of cone inputs to
% cells of the parvocellular pathway investigated with cone-isolating
% gratings", JOSA, 29, 2.
function BPI = STFbandpassIndex(theSpatialFrequencySupport, theSTF)

    Rmax = max(theSTF(:));
    [~,idx] = min(theSpatialFrequencySupport);
    Ro = theSTF(idx);
    BPI = Ro/Rmax;

end

