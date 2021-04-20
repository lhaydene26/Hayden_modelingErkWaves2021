% model parameters
alpha_1 = 12.6;
beta_1 = 0.35;
gamma_1 = 15.4;
alpha_2 = 0; % this is zero because we are outside the source region
alpha_3 = 3.9;
gamma_2 = 11.8;
alpha_4 = .5; 
gamma_3 = 0.14;
gamma_e = 0.14 + [0.0130,0.0143,0.0378,0.0804,0.1152,0.1401,0.1561,0.1642,0.1644]*gamma_1; % the matrix is the values that I takes along the leading edge each hour, taken manually from the standard simulation

% generate a grid of points for each nullcline in the leading edge (I = constant)
grid=0:.0001:1;
nullcline_1=zeros(numel(gamma_e),numel(grid));
nullcline_2=zeros(numel(gamma_e),numel(grid));
for i = 1:numel(gamma_e)
    for j=1:numel(grid)
        nullcline_1(i,j) = ((alpha_1*grid(j)^2)/(beta_1^2 + grid(j)^2))/(gamma_e(i) + (alpha_1*grid(j)^2)/(beta_1^2 + grid(j)^2));
        nullcline_2(i,j) = (gamma_2*grid(j) - alpha_2)/alpha_3;
    end
end

%  find the intersections of the nullclines
intersect_coords = cell(size(nullcline_1,1),1);
for i = 1:size(nullcline_1,1)
    [X0,Y0] = intersections(grid,nullcline_1(i,:),grid,nullcline_2(i,:),0);
    intersect_coords{i,1} = X0;
    intersect_coords{i,2} = Y0;
end


% generate the figure
fig = figure;
hold on;
n_plot1 = nullcline_1(1:2:end-1,:);
cmap = cool(size(n_plot1,1));
for i = 1:size(n_plot1,1)
    plot(grid,n_plot1(i,:),'color',cmap(i,:),'LineWidth',3);
end
plot(grid,nullcline_2(1,:),'k','LineWidth',3);
axis([0, 0.5, 0, 1]);
axis square
xlabel('[Growth Factor] (a.u.)') % x-axis label
ylabel('Active Erk (a.u.)') % y-axis label
legend('0h', '2h', '4h', '6h','Location','southeast');
legend boxoff
FigName = 'figures/leading_edge_nullclines';
standardizePlot(gcf,gca,FigName);
close(fig);