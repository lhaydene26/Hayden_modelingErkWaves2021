% This function loads the data generated by diffusible_inhibitor_model_series_xxx.m scripts for xx = 1:441,
% determines if stationary waves are formed, and calculates wave speed if normal waves are formed

clear data_store
data_store = cell(441,1);
data_store_whole = cell(441,1);
for i = 1:441
    if (exist(sprintf(strcat('diffusible_inhibitor_series_%d.txt'),i), 'file') == 2)
        data_store{i} = importdata(sprintf(strcat('diffusible_inhibitor_series_%d.txt'),i));
        data_store_whole{i} = importdata(sprintf(strcat('diffusible_inhibitor_series_whole_%d.txt'),i));
    end
end

dvary = linspace(0,10,21);
D_A = NaN(441,1);
D_I = NaN(441,1);
for i = 1:441
    D_A(i) = 566*dvary(mod(i-1,numel(dvary))+1); % diffusion constant for A
    D_I(i) = 566/10*3*dvary(ceil(i/numel(dvary))); % diffusion constant for I
end

tprint = 0.05;
tlength = 100;
tvals = 0:tprint:tlength;

xprint = 2;
xlength = 600;
xvals = 0:xprint:xlength;


% In lieu of a computational method, this section of code opens the last frame of each simulation and prompts
% the user to determine if stationary waves are formed

% load('workspaces\isNotStationary');
% isNotStationary = NaN(numel(data_store_pic),1);
% i = 1;
%  while (i <= numel(isNotStationary))
%     fig = figure;
%     imagesc(data_store_pic{i}, [0, 0.6]);
%     isNotStationary(i) = input('1 if no stationary waves, 0 if there are stationary waves');
%     close(fig)
%     if (isNotStationary(i) == 2)
%         i = i - 1;
%     end
%     i = i + 1;
% end
% isNotStationary = logical(isNotStationary);
% save('workspaces\isNotStationary','isNotStationary','-v7.3');


% calculate wave speed
wave_num = 1;
wave_distance = 600;
velocity_vals = NaN(21,21);
for i = 1:numel(data_store)
    if (and(numel(data_store{i}) > 0, isNotStationary(i)))
        [~, locs_l] = findpeaks(spline(tvals,data_store{i}(:,(wave_distance-20)/xprint+1),0:tprint/5:tlength),'MinPeakProminence',0.01);
        [~, locs_R] = findpeaks(spline(tvals,data_store{i}(:,(wave_distance-10)/xprint+1),0:tprint/5:tlength),'MinPeakProminence',0.01);
        [~, locs_end] = findpeaks(spline(tvals,data_store{i}(:,wave_distance/xprint+1),0:tprint/5:tlength),'MinPeakProminence',0.01);
        [~, locs_end2] = findpeaks(spline(tvals,data_store{i}(:,wave_distance/xprint+1),0:tprint/5:tlength),'MinPeakProminence',0.5);
        
        if (and(numel(locs_end) > 0,abs(numel(locs_end) - numel(locs_end2)) < 2))
            v1 = (10)/(locs_R(wave_num) - locs_l(wave_num))/tprint*5;
            v2 = (10)/(locs_end(wave_num) - locs_R(wave_num))/tprint*5;
            velocity_vals(i) = (v1+v2)/2;
        end
    end
    i
end


% make velocity/period plot
% fig = figure;
% imagesc(velocity_vals)
% ylabel('D_{activator}','FontSize',36);
% % ylim([0, 50]);
% xlabel('D_{inhibitor}','FontSize',36);
% colorbar
% xticks([0,5,10,15,20])
% xticklabels({'0','1','2','3','4'});
% yticks([0,5,10,15,20]);
% standardizePlot(gcf,gca,strcat('figures/DIvsDAvsV'));
% close(fig);


x = linspace(0,3,21);
x2 = [x(1)-(x(2)-x(1)),x];
[nr,nc] = size(velocity_vals);
fig = figure;
a = pcolor(0:3/20:3+3/20 ,0:0.5:10.5, [velocity_vals nan(nr,1); nan(1,nc+1)]);
shading flat;
ylabel('D_{activator}','FontSize',36);
xlabel('D_{inhibitor}','FontSize',36);
colorbar
xticks([0,0.5,1,1.5,2])
yticks([0,2,4,6,8,10]);
xlim([0,2]);
c = colorbar;
c.Label.String = 'Wave peak speed (um/h)';
standardizePlot_colorbar(gcf,gca,strcat('DIvsDAvsSpeed'));
close(fig);