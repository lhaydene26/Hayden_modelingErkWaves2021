% This is a representative function which runs a 2D simulation of the constant
% deposition model with diffusible inhibitor and varies the diffusion constants
% of the activator and inhibitor independently
%
% Scripts with the following variable: number = 1-441 were created and renamed
% appropriately then run on a computing cluster.
%
% Each of these scripts will generate a text file which is then analyzed
% with the script: analyze_diffusible_inhibitor_series.m

number = 1;
dvary = linspace(0,10,21);

D_A = 566*dvary(mod(number-1,numel(dvary))+1); % diffusion constant for the activator
D_I = 566/10*3*dvary(ceil(number/numel(dvary))); % diffusion constant for the inhibitor

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
D_coeff = [0 1 0; 1 -4 1; 0 1 0]; % matrix for diffusion

a2_source = 0.112*5; % value for a2 in the activator source region, chosen higher than normal to start the first wave sooner

% Space parameters
dx = 2; % space resolution for calculations
scaleSize = 500/dx; % radius of the physical scale object
N_pts = scaleSize*2 + 50/dx; % number of space units for the grid

% initial conditions
E = zeros(N_pts, N_pts); % initial values for E
A = zeros(N_pts, N_pts); % initial values for G
I = zeros(N_pts, N_pts); % initial values for I

% defining regions of space
xc = floor(N_pts/2); % center of the grid
yc = floor(N_pts/2); % center of the grid
[XX, YY] = meshgrid(1:size(E,2), 1:size(E,1)); % grid for the x and y coords
regionCX = xc-200/dx; % X center of the activator source region
regionCY = yc;% Y center of the activator source region
scaleRegion = sqrt((XX-xc).^2 + (YY-yc).^2) < scaleSize;
dRegion = 40/dx; % radius of the G source region
sourceRegion = sqrt((XX-regionCX).^2 + (YY-regionCY).^2) < dRegion;

% time parameters
timelength = 100; %total time
dt = 0.00002;
tprint = 0.05; % time interval to write to file
t_a2source = timelength; % time at which the source stops producing activator

xlength = 301;
resvals = NaN(timelength/tprint+1,xlength);

i = 1;
for time = 0:dt:timelength
    if (mod(time,tprint) == 0)        
        for j = 1:xlength
            resvals(i,j) = E(regionCY,regionCX + 1*(j-1));
        end
        i = i + 1;
        
    end
    
    [VE, VA, VI] = funct_potential(E,A,I,params);
  
    if (time < t_a2source)
        VA(sourceRegion) = VA(sourceRegion) + a2_source;
    end
    
    dE = (VE).*dt;
    dA = (VA+D_A.*funct_finiteDer_2D(A,D_coeff)./dx.^2).*dt;
    dI = (VI+D_I.*funct_finiteDer_2D(I,D_coeff)./dx.^2).*dt;
    
    E = E + dE;
    A = A + dA;
    I = I + dI;
    
    %boundary conditions
    E(~scaleRegion) = 0;

end

file_name = sprintf('diffusible_inhibitor_series_%d.txt',number);
fileID = fopen(file_name,'w');
for i = 1:size(resvals,1)
    fprintf(fileID,'%g\t',resvals(i,:));
    fprintf(fileID,'\n');
end
fclose(fileID);

resvals = E;
file_name = sprintf('diffusible_inhibitor_series_whole_%d.txt',number);
fileID = fopen(file_name,'w');
for i = 1:size(resvals,1)
    fprintf(fileID,'%g\t',resvals(i,:));
    fprintf(fileID,'\n');
end
fclose(fileID);