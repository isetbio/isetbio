%% Matlab convolution (conv2) and the photocurrent response
%
% These notes explain how we set up the conv2 in
%
%    @osLinear.computeCurrent
%
% In computeCurrent we use conv2(h,A) to calculate the photocurrent.  The
% columns of A are the absorptions, and h is the convolution kernel for a
% particular type of cone (L,M or S).
%
% The lesson illustrated here is that the impulse response function, h, in
% the conv2(h,A) should be extrapolated to times that match the full set of
% times of the abosrptions in A.
%
% Also, we trunctate the response at the end of the time period in which
% photon absorptions are specified.
%
% The reason for both of these choices are explained below.
%
% BW, ISETBIO Team, 2016

%% Suppose this is the convolution kernel

% We make an asymmetric kernel because it is easier to understand

tSamples = 1:100;
h = exp(-1*tSamples/10);
h = h(:);
plot(h);

%% Now, suppose we have a 2D data set (think photocurrent) here.

% Each column is a pixel, imagine the absorptions over time

% All zeros
A = zeros(100,3);

% Except one absorption
A(3,:) = 1;

%% When we convolve we get this

r = conv2(h,A);

fprintf('Length of r %d\n',length(r))

% A 1 in the first slot creates the exponential impulse response function.
% It last a long time and it doesn't start until a few steps in because,
% well, the absorptions start a few steps in
vcNewGraphWin; plot(r(:,1))

%% A value in the 3rd and 30th slot produces two of them

A = zeros(100,3);
A(3,:) = 1;
A(30,:) = 1;
r = conv2(h,A);
vcNewGraphWin; plot(r(:,1))

%%  Suppose we shorten the duration of the impulse response 

% If the IR is shorter than the duration of A ....
tSamples = 1:20;
h = exp(-1*tSamples/10);
h = h(:);

A = zeros(100,3);
A(3,:) = 1;
A(30,:) = 1;
r = conv2(h,A);

% The output is truncated.
vcNewGraphWin; plot(r(:,1))


%% But if the impulse response is longer than the stimulus

% We seem to be OK
tSamples = 1:100;
h = exp(-1*tSamples/10);
h = h(:);

A = zeros(30,3);
A(3,:) = 1;
A(30,:) = 1;
r = conv2(h,A);

% The output is truncated.
vcNewGraphWin; plot(r(:,1))

%% We decided to truncate the response to the length of the absorptions

% conv2 produces currents that carry on past the end of the stimulus time
% (absorptions).  
%
% BUT when a user creates a stimulus we think the user should bes
% specifying the entire time of interest, say by appending a uniform field.
% The values calculated by conv2() beyond the end of the stimulus assume
% that there are 0 photon absorptions in the stimulus. We don't think this
% zero-padding is expected by most (naive) users.

% We seem to be OK
tSamples = 1:100;
h = exp(-1*tSamples/10);
h = h(:);

A = zeros(30,3);
A(3,:) = 1;
A(30,:) = 1;
r = conv2(h,A);

r = r(1:100,:);

% The output is truncated.
vcNewGraphWin; plot(r(:,1))
