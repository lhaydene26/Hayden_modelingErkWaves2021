clear dataE dataA dataI
for i = 0:333
    if (exist(sprintf('DCC_text_files/dep_maxRexc_waveDist_E_%d.txt',i), 'file') == 2)
        dataE{i+1} = importdata(sprintf('DCC_text_files/dep_maxRexc_waveDist_E_%d.txt',i));
        dataA{i+1} = importdata(sprintf('DCC_text_files/dep_maxRexc_waveDist_A_%d.txt',i));
        dataI{i+1} = importdata(sprintf('DCC_text_files/dep_maxRexc_waveDist_I_%d.txt',i));
    end
end
save('workspaces/dataE.mat','dataE','-v7.3');
save('workspaces/dataA.mat','dataA','-v7.3');
save('workspaces/dataI.mat','dataI','-v7.3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tprint = 0.1;
tlength = 50;
tvals = 0:tprint:tlength;

dx = 1;
dstep = 1;
xmax = 200;
xvals = 0:dx*dstep:xmax*dx*dstep;

a1 = 12.6;
b1 = 0.35;
g1 = 15.4;
ge = 0.14;
a3 = 3.9;
g2 = 11.8;
a4 = 0.5;
g3 = 0.14;

A_vals = [0:0.02:0.2, 0.21:0.01:0.30, 0.301:0.001:0.361, 0.36101:0.00001:0.362, 0.3621:0.0001:0.363, 0.364:0.001:0.37, 0.39:0.02:0.5]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pwv = 12.1; % planar wave velocity of the first wave calculated from the planar wave velocity script
pwv = 12.7; % rough planar wave velocity calculated from the wave velocity in this experiment
D = 566; % diffusion constant
R_crit = round(D/pwv);

% go through all time points for the first 20 hours, in 20 hours, the decision to propagate or die off will have been decided
maxRexc = zeros(numel(dataE),3);
for i = 1:numel(dataE)
    tempE = dataE{i};
    tempA = dataA{i};
    tempI = dataI{i};
    Rexc = zeros(351,1);
    for j = 1:351
        dist_fE = tempE(j,:);
        dist_fA = tempA(j,:);
        dist_fI = tempI(j,:);
   
        for k = 1:numel(dist_fE)
            % determine if the point is excitable
            dt = 0.01;
            tlength = 20;
            EAI = zeros(tlength/dt+1,3);
            EAI(1,1) = dist_fE(k);
            EAI(1,2) = dist_fA(k);
            EAI(1,3) = dist_fI(k);
            times = 0:dt:tlength;
            for ii = 2:numel(times)
                EAI(ii,1) = EAI(ii-1,1) + dt*((a1*EAI(ii-1,2)^2)/(b1^2 + EAI(ii-1,2)^2)*(1-EAI(ii-1,1)) - EAI(ii-1,1)*(g1*EAI(ii-1,3) + ge));
                EAI(ii,2) = EAI(ii-1,2) + dt*(a3*EAI(ii-1,1) - g2*EAI(ii-1,2));
                EAI(ii,3) = EAI(ii-1,3) + dt*(g3*(a4*EAI(ii-1,1) - EAI(ii-1,3)));
            end
            if (or(max(EAI(:,1)) >= 0.2,max(tempE(1:j,k))>0.6))
                Rexc(j) = Rexc(j) + 1;
            else
                break;
            end
        end
    end
    maxRexc(i,2) = min(max(Rexc),100);
    
    % now calculate the farthest distance that the wave traveled
    thresh = 0.1; % threshold at or above which the wave is considered to be present at a location
    j = 101; % start at 100um away from the source. assumption that if the wave dies out, it won't reach that far before doing so. assumption based on looking at boundary near "almost" propagating waves
    while (j > 0)
        [pks, locs] = findpeaks(dataE{i}(:,j));
        if (and(numel(pks) > 0,max(pks) >= thresh))
            maxRexc(i,3) = (j-1)*dx;
            break;
        else
            j = j - 1;
        end
    end
    maxRexc(i,1) = A_vals(i);
    i
end

fig = figure;
times = [0, 5, 10, 15];
cmap = cool(numel(times));
hold on;
for i = 1:numel(times)
    plot(xvals,dataA{80}(times(i)/tprint+1,:)/max(dataA{80}(:)),'color',cmap(i,:),'LineWidth',3)
end
axis([0, 100, 0, 0.2])
xlabel('Radial distance (um)','FontSize',20)
ylabel('[Activator] (a.u.)','FontSize',20);
legend({'0h','5','10','15'});
legend boxoff
standardizePlot(gcf,gca,'4E_schematic');
close(fig);


fig = figure;
cmap = lines(size(maxRexc,1));
plot(maxRexc(:,2), maxRexc(:,3),'o','LineWidth',3, 'MarkerSize',20);
hold on;
plot([R_crit, R_crit], [0, 100],'k--','LineWidth',3);
xlabel('Maximum excited radius (um)','FontSize',36) % x-axis label
set(gca,'YTick',[0, 20, 40, 60, 80, 100])
set(gca,'YTickLabel',{'0','20','40','60','80','>100'})
set(gca,'XTick',[0, 20, 40, 80])
set(gca,'XTickLabel',{'0','20','40','stable wave'})
ylabel('Distance wave travels (um)','FontSize',36) % y-axis label
legend off
standardizePlot_wideAxes(gcf,gca,strcat('maxRexcited_vs_waveDistance'));
close(fig);
