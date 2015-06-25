% t_outersegment_short
%A very brief tutorial introducing the outersegment object
% This will be updated June 26
% James Golden

clear

% load sensor
load('sensor_ex.mat');

noise_flag = 1;
os1 = outersegment(noise_flag);

os1 = outersegmentCompute(os1, sensor, 'linearfilter');

os1

os2 = outersegment(noise_flag);

os2 = outersegmentCompute(os1, sensor, 'nonlinear1');

os2