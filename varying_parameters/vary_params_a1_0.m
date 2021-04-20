% This is a representative function which varies the value of alpha_1. 
% Scripts with the following variable: number = 0-100 were created and renamed
% appropriately then run on a computing cluster.
%
% Each of these scripts will generate a text file which is then analyzed
% with the script "analyze_xx.m" where xx is the parameter shorthand (alpha_1 = a1, etc.)

number = 0;
vary_vals = linspace(0, 3, 101);

% model parameters
params = [];
params.alpha_1 = 12.6*vary_vals(number+1);
params.beta_1 = 0.35;
params.gamma_1 = 15.4;
params.alpha_2 = 0;
params.alpha_3 = 3.9;
params.gamma_2 = 11.8;
params.alpha_4 = .5; 
params.gamma_3 = 0.14;
params.gamma_e = 0.14;
a2_source = 0.112; % value for a2 in the activator source region
D = 566; % diffusion constant

dx = 1; % space resolution for calculations
D_coeff = [0 1 0; 1 -4 1; 0 1 0]; % matrix for diffusion

% Space parameters
N_pts_y = 10/dx; % number of space units for the grid x
N_pts_x = 1000/dx;

% initial conditions
E = zeros(N_pts_y, N_pts_x); % initial values for E
A = zeros(N_pts_y, N_pts_x); % initial values for G
I = zeros(N_pts_y, N_pts_x); % initial values for I

% defining regions of space, using a long strip as our domain to generate planar waves
[XX, YY] = meshgrid(1:size(E,2), 1:size(E,1)); % grid for the x and y coords
regionCX = 1; % X center of the activator source region
dRegion = 50/dx; % radius of the activator source region
sourceRegion = XX <= dRegion;

% time parameters
timelength = 200; %total time
dt = 0.0001; % time resolution for calcuations
tprint = 0.1;

resvals = zeros(timelength/tprint+1,201);
i = 1;
for time = 0:dt:timelength
    if (mod(time,tprint) == 0)
        for j = 1:201 % only captuure the first 200 um unless varying diffusion which increases wave width substantially
            resvals(i,j) = E(5,regionCX + 1*(j-1));
        end
        i = i + 1;
    end
  
    [VE, VA, VI] = funct_potential(E,A,I,params);
    

    VA(sourceRegion) = VA(sourceRegion) + a2_source;

    dE = (VE).*dt;
    dA = (VA+D.*funct_finiteDer_2D_symmetric(A,D_coeff)./dx.^2).*dt;
    dI = (VI).*dt;
    
    E = E + dE;
    A = A + dA;
    I = I + dI;
    
end

file_name = sprintf('vary_params_a1_%d.txt',number);
fileID = fopen(file_name,'w');
for i = 1:size(resvals,1)
    fprintf(fileID,'%g\t',resvals(i,:));
    fprintf(fileID,'\n');
end
fclose(fileID);