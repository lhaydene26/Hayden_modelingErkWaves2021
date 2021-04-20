% This code is a 2D simulation of the oscillatory deposition model

% model parameters
params = [];
params.alpha_1 = 12.6;
params.beta_1 = 0.35;
params.gamma_1 = 15.4;
params.alpha_2 = 0;
params.alpha_3 = 3.9;
params.gamma_2 = 11.8;
params.alpha_4 = .5; 
params.gamma_3 = 0.14;
params.gamma_e = 0.14;
D = 566; % diffusion constant
D_coeff = [0 1 0; 1 -4 1; 0 1 0]; % matrix for diffusion

a = 2.8;
s = 1;
T = 24.3;
ph = T/2;

% Space parameters
dx = 5; % space resolution for calculations
scaleSize = 520/dx; % radius of the physical scale object
N_pts = scaleSize*2 + 50/dx; % number of space units for the grid

% initial conditions
E = zeros(N_pts, N_pts); % initial values for E
A = zeros(N_pts, N_pts); % initial values for G
I = zeros(N_pts, N_pts); % initial values for I

% defining regions of space
xc = floor(N_pts/2); % center of the grid
yc = floor(N_pts/2); % center of the grid
[XX, YY] = meshgrid(1:size(E,2), 1:size(E,1)); % grid for the x and y coords
regionCX = xc-N_pts/10; % X center of the G source region
regionCY = yc-N_pts/10; % Y center of the G source region
scaleRegion = sqrt((XX-xc).^2 + (YY-yc).^2) < scaleSize;
dRegion = 40/dx; % radius of the G source region
sourceRegion = sqrt((XX-regionCX).^2 + (YY-regionCY).^2) < dRegion;

% time parameters
timelength = 100; %total time
dt = 0.005; % time resolution for calcuations
tprint = 0.1; % time interval to write to file
t_a2source = timelength; % time at which the source stops producing activator

% video settings
vid = VideoWriter('oscillatory_deposition_model','MPEG-4');
vid.FrameRate = 12;
open(vid);

for time = 0:dt:timelength
    if (mod(time,tprint) == 0)
        
        fig = figure(1);
        finalBox = [360, 360]; % for simulations
        imgInit = E;
        finalView = imref2d(finalBox,[size(imgInit,2)/2-finalBox(2)/2 size(imgInit,2)/2+finalBox(2)/2],[size(imgInit,1)/2-finalBox(1)/2 size(imgInit,1)/2+finalBox(1)/2]);
        imgFinal = imwarp(imgInit,affine2d([1 0 0; 0 1 0; 0 0 1]),'OutputView',finalView);
        scaleRegion_Final = imwarp(scaleRegion,affine2d([1 0 0; 0 1 0; 0 0 1]),'OutputView',finalView);
        img = imagesc(imgFinal, [0, 0.6]); hold on;
        img.AlphaData = scaleRegion_Final;
        
        axis equal
        axis tight
        set(gcf, 'InvertHardCopy', 'off');
        mygcf=gcf;
        mygcf.Units = 'centimeters';
        mygcf.Position = [mygcf.Position(1) mygcf.Position(2)-5 16 16];
        ax=gca;
        set(gca,'Color','k')
        set(gca,'XTick',[])
        set(gca,'YTick',[])
        ax.Position=[0 0 1 1];
        text(0.8, 0.92, strcat(num2str(time, '%.1f'),' hr'), 'Units', 'Normalized', 'Color', 'w', 'FontSize', 40, 'FontWeight', 'bold','HorizontalAlignment','right');
        drawnow;
        
        f = getframe;
        cdata2 = imresize(f.cdata,1);
        f.cdata = cdata2;
        writeVideo(vid,f);
        close(fig);
    end
    
    [VE, VA, VI] = funct_potential(E,A,I,params);

    if (time < t_a2source)
        VA(sourceRegion) = VA(sourceRegion) + a/(sqrt(2*pi())*s)*exp(-(mod(time,T)-ph)^2/(2*s^2));
    end
    
    dE = (VE).*dt;
    dA = (VA+D.*funct_finiteDer_2D(A,D_coeff)./dx.^2).*dt;
    dI = (VI).*dt;
    
    E = E + dE;
    A = A + dA;
    I = I + dI;
    
    %boundary conditions for E and I
    E(~scaleRegion) = 0;
    I(~scaleRegion) = 0;

end

close(vid)