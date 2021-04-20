% This is a representative function which runs a 2D simulation of the constant
% deposition model with a variable sized activator pool in the source region and no additional production
%
% Scripts with the following variable: number = 1-333 were created and renamed
% appropriately then run on a computing cluster.
%
% Each of these scripts will generate a text file which is then analyzed
% with the script: analyze_maxRexc_vs_wave_distance.m

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
dx = 1; % space resolution for calculations
scaleSize = 400/dx; % radius of the physical scale object
N_pts = scaleSize*2 + 50/dx; % number of space units for the grid

% initial conditions
E = zeros(N_pts, N_pts); % initial values for E
A = zeros(N_pts, N_pts); % initial values for A
I = zeros(N_pts, N_pts); % initial values for I

% defining regions of space
xc = floor(N_pts/2); % center of the grid
yc = floor(N_pts/2); % center of the grid
[XX, YY] = meshgrid(1:size(E,2), 1:size(E,1)); % grid for the x and y coords
regionCX = floor(xc-200/dx); % X center of the source region
regionCY = yc; % Y center of the source region
scaleRegion = sqrt((XX-xc).^2 + (YY-yc).^2) < scaleSize;
dRegion = 10/dx; % radius of the source region
sourceRegion = sqrt((XX-regionCX).^2 + (YY-regionCY).^2) < dRegion;
A_vals = [0:0.02:0.2, 0.21:0.01:0.30, 0.301:0.001:0.361, 0.36101:0.00001:0.3613, 0.361301:0.0000005:0.3614, 0.3614:0.00001:0.3612, 0.3621:0.0001:0.363, 0.364:0.001:0.37, 0.39:0.02:0.5]'; % a manually determiend range which spans the response range of wave behaviors
A_init = A_vals(number+1);
A(sourceRegion) = A_init;

% time parameters
timelength = 50; %total time
dt = 0.0002; % time resolution for calcuations
tprint = 0.1;

resvalsE = zeros(timelength/tprint+1,201);
resvalsA = zeros(timelength/tprint+1,201);
resvalsI = zeros(timelength/tprint+1,201);
i = 1;
for time = 0:dt:timelength
    if (mod(time,tprint) == 0)
        for j = 1:201
            resvalsE(i,j) = E(regionCY,regionCX + 1*(j-1));
            resvalsA(i,j) = A(regionCY,regionCX + 1*(j-1));
            resvalsI(i,j) = I(regionCY,regionCX + 1*(j-1));
        end
        i = i + 1;
    end
    
    [VE, VA, VI] = funct_potential(E,A,I,params);
    
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

file_name = sprintf('dep_maxRexc_waveDist_E_%d.txt',number);
fileID = fopen(file_name,'w');
for i = 1:size(resvalsE,1)
    fprintf(fileID,'%g\t',resvalsE(i,:));
    fprintf(fileID,'\n');
end
fclose(fileID);

file_name = sprintf('dep_maxRexc_waveDist_A_%d.txt',number);
fileID = fopen(file_name,'w');
for i = 1:size(resvalsA,1)
    fprintf(fileID,'%g\t',resvalsA(i,:));
    fprintf(fileID,'\n');
end
fclose(fileID);

file_name = sprintf('dep_maxRexc_waveDist_I_%d.txt',number);
fileID = fopen(file_name,'w');
for i = 1:size(resvalsI,1)
    fprintf(fileID,'%g\t',resvalsI(i,:));
    fprintf(fileID,'\n');
end
fclose(fileID);