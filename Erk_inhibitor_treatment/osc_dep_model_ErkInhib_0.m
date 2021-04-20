% This is a representative function which runs a 2D simulation of the oscillatory
% deposition model and varies the time at which an Erk inhibitor is added and removed.
%
% Scripts with the following variable: number = 0-100 were created and renamed
% appropriately then run on a computing cluster.
%
% Each of these scripts will generate a text file which is then analyzed
% with the script: analyze_osc_dep_model_ErkInhib.m

number = 0;

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
scaleSize = 400/dx; % radius of the physical scale object
N_pts = scaleSize*2 + 50/dx; % number of space units for the grid

% initial conditions
E = zeros(N_pts, N_pts); % initial values for E
A = zeros(N_pts, N_pts); % initial values for G
I = zeros(N_pts, N_pts); % initial values for I

% defining regions of space
xc = floor(N_pts/2); % center of the grid
yc = floor(N_pts/2); % center of the grid
[XX, YY] = meshgrid(1:size(E,2), 1:size(E,1)); % grid for the x and y coords
regionCX = xc; % X center of the G source region
regionCY = yc; % Y center of the G source region
scaleRegion = sqrt((XX-xc).^2 + (YY-yc).^2) < scaleSize;
dRegion = 50/dx; % radius of the G source region
sourceRegion = sqrt((XX-regionCX).^2 + (YY-regionCY).^2) < dRegion;

% time parameters
timelength = 150; %total time
dt = 0.005; % time resolution for calcuations
tprint = 0.1;

t_washout = number; % time at which the drug treatment starts
washout_duration = 10; % how long the treatment lasts

resvals = zeros(timelength/tprint+1,2);
i = 1;
for time = 0:dt:timelength
    if (mod(time,tprint) == 0)
        resvals(i,1) = E(regionCY,regionCX);
        resvals(i,2) = E(regionCY,regionCX + 200/dx);
        i = i + 1;
    end
    
    [VE, VA, VI] = funct_potential(E,A,I,params);
    VA(sourceRegion) = VA(sourceRegion) + a/(sqrt(2*pi())*s)*exp(-(mod(time,T)-ph)^2/(2*s^2));
    
    params.alpha_1 = 12.6;
    if (time >= t_washout)
        if (time < t_washout+washout_duration)
            params.alpha_1 = 0;
        end
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

file_name = sprintf('osc_dep_ErkInhib_%d.txt',number);
fileID = fopen(file_name,'w');
for i = 1:size(resvals,1)
    fprintf(fileID,'%g\t',resvals(i,:));
    fprintf(fileID,'\n');
end
fclose(fileID);