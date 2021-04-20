% This code is a 2D simulation of the constant deposition model and
% additionally plots the activator concentration over space at different
% points in time


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
a2_source = 0.112; % value for a2 in the G source region
D = 566; % diffusion constant
dx = 5; % space resolution for calculations
D_coeff = [0 1 0; 1 -4 1; 0 1 0]; % matrix for diffusion

% Space parameters
N_pts = 1600/dx; % number of space units for the grid


% initial conditions
E = zeros(N_pts, N_pts); % initial values for E
A = zeros(N_pts, N_pts); % initial values for A
I = zeros(N_pts, N_pts); % initial values for I

% defining regions of space
xc = floor(N_pts/2); % center of the grid
yc = floor(N_pts/2); % center of the grid
[XX, YY] = meshgrid(1:size(E,2), 1:size(E,1)); % grid for the x and y coords
regionCX = floor(xc-300/dx); % X center of the activator source region
regionCY = yc; % Y center of the activator source region
dRegion = 50/dx; % radius of the activator source region
sourceRegion = sqrt((XX-regionCX).^2 + (YY-regionCY).^2) < dRegion;

scaleSize = 750/dx; % radius of the physical scale object
scaleRegion = sqrt((XX-xc).^2 + (YY-yc).^2) <= scaleSize;

% time parameters
timelength = 200; %total time
dt = 0.001; % time resolution for calcuations
tprint = 0.01;
tsource = 10.5;

resvals = zeros(timelength/tprint+1,101);
i = 1;
for time = 0:dt:timelength
    if (mod(time,tprint) == 0)
        for j = 1:101
            resvals(i,j) = A(regionCY,regionCX + 10/dx*(j-1));
        end
        i = i + 1;
    end
    
    [VE, VA, VI] = funct_potential(E,A,I,params);
    if (time < tsource)
        VA(sourceRegion) = VA(sourceRegion) + a2_source;
    end
    
    dE = (VE).*dt;
    dA = (VA+D.*funct_finiteDer_2D(A,D_coeff)./dx.^2).*dt;
    dI = (VI).*dt;
    
    E = E + dE;
    A = A + dA;
    I = I + dI;
    
    %boundary conditions
    E(~scaleRegion) = 0;
    I(~scaleRegion) = 0;
end

xlength = 1000;
tvals = 0:tprint:timelength;
xvals = 0:10:xlength;

tcmap = cool(5);

figure;
before = findall(gca);
fnplt(spapi(5,xvals,resvals(1/dt*6+1,:)),2);
added = setdiff(findall(gca), before);
set(added, 'Color', tcmap(1,:))
hold on
before = findall(gca);
fnplt(spapi(5,xvals,resvals(1/dt*12+1,:)),2);
added = setdiff(findall(gca), before);
set(added, 'Color', tcmap(2,:))
hold on
before = findall(gca);
fnplt(spapi(5,xvals,resvals(1/dt*24+1,:)),2);
added = setdiff(findall(gca), before);
set(added, 'Color', tcmap(3,:))
hold on
before = findall(gca);
fnplt(spapi(5,xvals,resvals(1/dt*48+1,:)),2);
added = setdiff(findall(gca), before);
set(added, 'Color', tcmap(4,:))
hold on
before = findall(gca);
fnplt(spapi(5,xvals,resvals(1/dt*72+1,:)),2);
added = setdiff(findall(gca), before);
set(added, 'Color', tcmap(5,:))
axis([0 1000 0 0.4])
set(gca,'FontSize',14)
set(gca,'LineWidth',2)
xlabel('Radial distance (um)','FontSize',20) % x-axis label
ylabel('[Activator] (a.u.)','FontSize',20) % y-axis label
legend('0.25 days','0.5','1','2','3','Location','NorthEast');
legend boxoff
FigName = strcat('wave_traces');
standardizePlot(gcf,gca,FigName);




