function rgcMovieWave(obj)
% 
% Generates a movie of the spikes in a mosaic of an rgcLayer object.

%           rgcMovieWave(obj)
% 
% This is an attempt at a new type of visualization of the spiking activity
% of an rgc mosaic. The spikes of each cell are plotted over time so that
% they appear to be moving out of the screen. The spiking information
% conveyed via lateral connections between cells is also shown for each
% spike as it travels to nearby cells.
% 
% Inputs: an rgc object.
% 
% Outputs: a movie showing the spiking of rgc cells over time.
% 
% 2/2016 JRG (c) isetbio team

%%
% Initialize figure
% vcNewGraphWin([],'upperleftbig');
figure; set(gcf,'position',[160 60 1070 740]);
hold on;

% TODO: Loop over all mosaics
nX = 0; nY = 0;
% for cellTypeInd = 1:5
cellTypeInd = 1;

[nX,nY,~] = size(obj.mosaic{cellTypeInd}.linearResponse);
% nX = nX + nXi;
% nY = nY = nYi;

% Set expansion factor of spike plots 
zfac = .002;

% Set frame subsampling number
frameskip= 10;
for t = 1:frameskip:5760
    hold on;
    cla
    
    % Loop through each cell and plot spikes over time 
    for xc = 1:nX
        for yc = 1:nY
            
            % Get the appropriate spike data
            spPlot=obj.mosaic{cellTypeInd}.spikeResponse{xc,yc,1,2}(t:t+1000);
            % spPlot=(median(horzcat(obj.mosaic{3}.spikeResponse{xc,yc,:,2})'));
            
            % Get the time values
            t1 = flipud((1:length(spPlot))');
            
            % xv sets x position of cell in spatial array
            xv = zfac*(xc-ceil(nX/2))*t1+(xc-ceil(nX/2))+zeros(length(spPlot),1);
            % yv sets time position using t
            yv = t/10000+1/10000*(1:length(spPlot));
            % zv sets y position of cell in array and scales size of spike
            % waveform
            spScale = 1;
            zv = 10*(yc-ceil(nY/2))+10*(yc-ceil(nY/2))*zfac*t1+spScale*spPlot;
            
            % Plot the waveform for this cell
            h1=plot3(xv, yv, zv,'linewidth',2);
            colorval = get(h1,'color');
            
            % Plot the spike activity transmitted by lateral connections
            % Find spikes in this temporal window
            % Get spike times for this cell
            spTimes = 100*obj.mosaic{cellTypeInd}.spikeResponse{xc,yc,1,1};
            % Check if any fall in the temporal window
            spFind = find(spTimes >= t & spTimes < t+frameskip);
            
            % If there are spikes, plot them on lateral connections
            if length(spFind)>0
                % Check which other cells are connected to the cell of
                % interest
                for xc2 = 1:nX
                    for yc2 = 1:nY
                        % If the coupling weight is > 0 and there is a
                        % spike, plot this
                        if abs(obj.mosaic{cellTypeInd}.couplingMatrix{xc,yc}(xc2,yc2))>0
                            % Calculate spatial positions in array of
                            % starting cell and target cell
                            x0 = xc-ceil(nX/2); xf = xc2-ceil(nX/2);
                            z0 = yc-ceil(nY/2); zf = yc2-ceil(nY/2);
                            hold on;
                            
                            % plot the line representing the lateral
                            % connection
                            line(1+zfac*1000*[x0 xf ],[t  t]/10000,1+zfac*1000*10*[z0 zf],'color',colorval,'linewidth',4);
                            
                        end%if
                    end%yc2
                end%xc2
                
            end%if length
            
        end%yc
    end%xc
    
    % Label axes
    xlabel('x position'); ylabel('time (sec)'); zlabel('y position');
    set(gca,'fontsize',14);
    
    % Set view angle
    % view(-19,18);
    view(-2,2)
    
    % Shift axis
    axis([-8 8 t/10000 (t+1000)/10000 -40 60]);
    grid on;
    drawnow;
    % pause(.01);
    % hold off;
    % clf
    
end

% Alternate method, plot whole waveform and shift axis
% for t = 1:10:6250
%     
% %     axis([-3 3 0+t 500+t -30 30])
%     axis([-4 4 0+t 1000+t -50 50])
%     pause(.005);
% end