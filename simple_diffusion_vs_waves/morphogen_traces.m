% This code is a 2D simulation of the constant deposition model with all parameters except a2 and D set to 0.
% This is analogous to a simple diffusion model where activator is produced and diffuses across the domain
% This code generates the data and plots the activator concentration over space at different
% points in time

% model parameters
a2_source = 0.112; % value for a2 in the G source region

D = 566; % diffusion constant
dx = 10; % space resolution for calculations
D_coeff = [0 1 0; 1 -4 1; 0 1 0]; % matrix for diffusion

% Space parameters
N_pts = 4000/dx; % number of space units for the grid

% initial conditions
A = zeros(N_pts, N_pts); % initial values for the activator

% defining regions of space
xc = floor(N_pts/2); % center of the grid
yc = floor(N_pts/2); % center of the grid
[XX, YY] = meshgrid(1:size(A,2), 1:size(A,1)); % grid for the x and y coords
regionCX = xc; % X center of the activator source region
regionCY = yc; % Y center of the activator source region
dRegion = 50/dx; % radius of the activator source region
sourceRegion = sqrt((XX-regionCX).^2 + (YY-regionCY).^2) < dRegion;

scaleSize = 2000/dx; % radius of the physical scale object
scaleRegion = sqrt((XX-xc).^2 + (YY-yc).^2) < scaleSize;

% time parameters
timelength = 24*50; %total time
dt = 0.01; % time resolution for calcuations
tprint = 6;

resvals = zeros(timelength/tprint+1,51);
i = 1;
for time = 0:dt:timelength
    if (mod(time,tprint) == 0)        
        for j = 1:51
            resvals(i,j) = A(regionCY,regionCX + 20/dx*(j-1));
        end
        i = i + 1;
    end
    
    VA = zeros(N_pts,N_pts);
    VA(sourceRegion) = VA(sourceRegion) + a2_source;
    
    dA = (VA+D.*funct_finiteDer_2D(A,D_coeff)./dx.^2).*dt;
    
    A = A + dA;
end


tcmap = cool(6);
fig = figure;
k = 1;
for i = [2,3,5,13,41,201]
    plot(0:20:1000,resvals(i,:),'Color',tcmap(k,:),'LineWidth',3)
    hold on
    k = k+1;
end
axis([0 1000 0 1])
set(gca,'FontSize',15)
set(gca,'LineWidth',2)
xlabel('Radial distance (um)','FontSize',20) % x-axis label
ylabel('[Morphogen] (a.u.)','FontSize',20) % y-axis label
legend('0.25 days','0.5','1','3','10','50','Location','NorthEast');
legend boxoff

FigName = strcat('morphogen_traces');
standardizePlot(gcf,gca,FigName);
close(fig);



