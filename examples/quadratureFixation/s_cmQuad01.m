%% Show the basic quadrature calculation
%
% These are all for the circular shift case, just checking the logic.
% Main conclusions:
%    When we use a whole harmonic, we get the expected result.  No
%    change as we circularly shift the pattern
%    When we use a Gabor, we get a different result.  There is a shift
%    that depends on some combination of the spatial frequency and the
%    size of the Gabor window
%
% Wandell

%%
ieInit
chdir(fullfile(isetbioRootPath,'local'));

%% Testing the basic idea

%  Make a signal and two harmonics in quadrature
x = (0:127)/128;
fquad = 5;
s = sin(2*pi*fquad*x);
c = cos(2*pi*fquad*x);

fsig = 5;
sig = 0.4*square(2*fsig*pi*x) + 0.5;

eBase = dot(sig,s)^2 + dot(sig,c)^2;

vcNewGraphWin; plot(x,sig,'k-',x,s,'r-',x,c,'b-');

%% Take the inner product of the signal with each harmonic. 

% Then compute the energy, also known as the amplitude at that
% frequency.

% Shift the signal and recompute
for ii=1:2:10
    eShift = dot(circshift(sig,ii),s)^2 + dot(circshift(sig,ii),c)^2;
    fprintf('Difference: %.6f\n',eBase - eShift)
end

%% Now do the same, but for a 2D image
img    = repmat(sig,[128,1]);
simg   = repmat(s,[128,1]);
cimg   = repmat(c,[128,1]);
vcNewGraphWin; imagesc(img); colormap(gray); axis image

eBase = dot(img(:),simg(:))^2 + dot(img(:),cimg(:))^2;
for ii=1:2:10
    thisIMG = circshift(img,ii,2);
    % imagesc(thisIMG); colormap(gray); axis image; pause(0.2);
    eShift = dot(thisIMG(:),simg(:))^2 + dot(thisIMG(:),cimg(:))^2;
    fprintf('Difference: %.6f\n',eBase - eShift)
end

%% Now apply by a Gaussian envelope, rather than using the harmonic

% Note:
%
% When fsig = fquad, and both are high, the Gaussian envelope has very
% little impact. But when fsig (frequency of the square wave) differs
% from fquad there is a big difference

% Big difference with a small envelope, and little difference with a big
% envelope, like the full harmonic above.
for spread = 16:16:64
    g = fspecial('gaussian',[128 128],spread);
    
    % vcNewGraphWin; imagesc(g); colormap(gray); axis image
    gsimg = g .* simg;
    gcimg = g .* cimg;
    
    eBase = dot(img(:),gsimg(:))^2 + dot(img(:),gcimg(:))^2;
    pError = zeros(spread+1,1);
    % Small shifts and there is a constant response
    fprintf('\n---- Spread %d\n',spread);
    for ii=0:spread
        thisIMG = circshift(img,ii,2);
        % imagesc(g.*thisIMG); colormap(gray); axis image; pause(0.2);
        eShift = dot(thisIMG(:),gsimg(:))^2 + dot(thisIMG(:),gcimg(:))^2;
        pError(ii+1) = 100*(eBase - eShift)/eBase;
        fprintf('Difference (percentage): %.6f (step %d)\n',pError(ii+1),ii);
    end
    vcNewGraphWin; 
    plot((0:spread),pError); title(sprintf('Spread %d',spread));
    line([0,spread],[5,5],'Color','k')
    line([0,spread],[-5,-5],'Color','k')
end

%% Make the simg and cimg this way?
%{
hparams = harmonicP; 
hparams.freq = 2; hparams.row = 128; hparams.col = 128;
simg = imageHarmonic(hparams);
simg = simg - 1;
vcNewGraphWin; mesh(simg);

hparams.ph = 0;
cimg = imageHarmonic(hparams);
cimg = cimg - 1;
vcNewGraphWin; mesh(cimg);
%}

%% END
