function ieMovie(movieMatrix)
% Show a movie of an (x,y,t) matrix with minimal fuss.

%% Show test movie

szMovie = size(movieMatrix);

vcNewGraphWin([],'upperleftbig'); 

for frame1 = 1:szMovie(3)
    imagesc(movieMatrix(:,:,frame1));
    caxis([0 1]);
    colormap gray; 
    drawnow;
end
close;