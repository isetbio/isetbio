%% Mask out specific cones
%
% Work with David to create a simple function that masks (zeros out?) 
% specific functions, perhaps defined as a spatial mask over the cone
% mosaic
%
%{
Hello Dr. Wandell, 

My name is David Buickians. I work in Dr. Moss's Lab, in the Ophthalmology Department.  I have been familiarizing myself with isetBio with isettools and working through the t_cones tutorials. Which have tremendously help in understanding the isetBio framework. 
 
My primary goal is to apply a mask on the cone mosaic such that those cones covered by the mask do not participate in the computation of the excitation or photon absorption.

I attempted to bootstrap already programmed parameters in cMosaic.m to accomplishing such a 2D mask. However, these parameters affect all the cones within the mosaic.
      1. eccVaryingConeAperture
      2. eccVaryingOuterSegmentLength
I then found a function within the scene folder labeled scenceCrop.m which I thought I could crop specific regions from an unchanged cMosaic which would mimic a masked cone mosaic. 

Questions: 
1. Is there a parameter in cMosiac that can accomplish a 2D mask?
      Is there a function that allows individual control of cone placement?
      Is there a function that allows individual control of cone absorption?
2. Should I create a cone mosaic which incorporates the mask into the mosaic and then use the cMaskMosaic or should I apply the mask after the scene has been applied to the cone mosaic matrix? 

Would it be possible to set a time to have a zoom meeting? 

Best, 
David Buickians
%}