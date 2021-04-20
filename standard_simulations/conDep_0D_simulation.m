% This script generates a 0D simulation of the constant deposition model.
% This is analogous to a simulation of a single point in space inside or outside the source region (as specified)

n = 2;
a1 = 12.6;
b1 = 0.35;
g1 = 15.4;
a2 = 0.112; % leave this uncommented to simulate inside the source region
% a2 = 0;   % leave this uncommented instead to simulate outside the source region
a3 = 3.9;
g2 = 11.8;
a4 = .5; 
g3 = 0.14;
ge = 0.14;

dt = 0.002;
tlength = 100;
E = NaN(tlength/dt+1,1);
E(1) = 0;
A = E;
I = E;
i = 2;
times = 0:dt:tlength;
for t = 1:numel(times)-1
    E(i) = E(i-1) + dt*(a1*A(i-1)^n/(b1^n+A(i-1)^n)*(1-E(i-1)) - g1*E(i-1)*I(i-1) - ge*E(i-1));
    A(i) = A(i-1) + dt*(a2 + a3*E(i-1) - g2*A(i-1));
    I(i) = I(i-1) + dt*(g3*(-I(i-1) + a4*E(i-1)));
    i = i + 1;
end 


fig = figure;
hold on;
plot(times,E,'color','#0A0AFF','LineWidth',3);
plot(times,A,'color','#009F11','LineWidth',3);
plot(times,I,'color','#E91300','LineWidth',3);
axis([0, 100, 0, 1]);
axis square
xlabel('Time (h)') % x-axis label
ylabel('Species activity (a.u.)') % y-axis label
legend('Erk', 'Activator', 'Inhibitor','Location','northeast');
legend boxoff
FigName = 'constant_deposition_model_0D';
% standardizePlot(gcf,gca,FigName);
% close(fig);