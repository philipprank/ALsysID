t1_r = squeeze(sys_id_rand(:,3,:));
t1_hz = squeeze(sys_id_rcdhz(:,3,:));
t1_tr = ones(36,100)*sys_sim.B(3);
v1 = (t0-t1)./t0*100;
v2 = (t0-t2)./t0*100;
v3 = abs(v2)-abs(v1);

%%
map = [0.0000 0.7451 1.0000; 0.0000 0.3176 0.6196;...
    0.2431 0.2667 0.2980; 0.6235 0.6000 0.5961];
f1 = figure;
%tiledlayout('flow')

%nexttile;
plot(randn(10,1),randn(10,1))
title('Hello')
xlabel('test')

% nexttile;
% plot(randn(40,1),randn(40,1))
% 
% nexttile;
% plot(randn(10,1),randn(10,1))

savefig(f1,"plot1.fig")


%%
surf(v3)
view(2)
colormap(map)
colorbar
clim([-60 20])
xlim([1 100])
ylim([1 35])
grid