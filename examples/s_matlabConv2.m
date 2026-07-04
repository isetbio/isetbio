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
% The impulse response function, h, in the conv2(h,A) is always returned to
% be at least 400 ms, and it is interpolated to match the times of the
% abosrptions in A.
%
% Also - we trunctate responses after the stimulus (photon absorptions) are
% specified.
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

A = zeros(200,3);
A(3,:) = 1;
A(30,:) = 1;
r = conv2(h,A);
vcNewGraphWin; plot(r(:,1))

%% A value in the 3rd and 30th slot produces two of them

A = zeros(60,3);
A(3,:) = 1;
A(30,:) = 1;
r = conv2(h,A);
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

%% We truncate the response to the length of the absorptions

%
% conv2 produces currents that carry on past the end of the stimulus time
% (absorptions).  
%
% BUT when a user creates a stimulus we think the user should sspecify the
% entire time of interest, say by appending a uniform field. The values
% calculated by conv2() beyond the end of the stimulus assume that there
% are 0 photon absorptions in the stimulus. We don't think this
% zero-padding is expected by most (naive) users.

% We seem to be OK
tSamples = 1:100;
h = exp(-1*tSamples/10);
h = h(:);

A = zeros(30,3);
A(3,:) = 1;
A(30,:) = 1;
r = conv2(h,A);

r = r(1:30,:);

% The output is truncated.
vcNewGraphWin; plot(r(:,1))

%%
