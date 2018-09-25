function supportDegs = supportInDegs(pixelsNum, pixelsDegs)
% Given the number of pixels and the spatial extent in degrees, compute
% the spatial support in degrees.
%
% 7/24/18  npc  Wrote it
%
    supportDegs = 1:pixelsNum;
    supportDegs = supportDegs - mean(supportDegs);
    supportDegs = 0.5*supportDegs / max(abs(supportDegs));
    supportDegs = supportDegs * pixelsDegs;
end