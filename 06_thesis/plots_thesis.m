%%================================%%
%%%%   Thesis Plots             %%%%
%%================================%%

%% Figure 3.1 - Probability Density Functions
layout_options;

f1 = figure;
f1.Units = 'centimeters';
f1.Position = [8 4 14 8];

mu_range = -2:0.1:10;
mu = 4;
sigma1 = 1;
sigma2 = 4;
obsv = 20;

y1 = normpdf(mu_range,mu,sigma1);
y2 = normpdf(mu_range,mu,sigma2);
plot(mu_range,y1,mu_range,y2, 'LineWidth',1.5)
xlabel('X')
ylabel('Probability Density $p(X)$')
ylim([0 0.42])
legend(sprintf('$\\sigma_{1} = %i$', sigma1),sprintf('$\\sigma_{2} = %i$', sigma2),...
    'Location','northeast')

exportgraphics(f1, 'fig_3_1.pdf', 'ContentType', 'vector');


%% Figure 3.2 - Likelihood Functions
data1 = normrnd(mu,sigma1,obsv,1);
data2 = normrnd(mu,sigma2,obsv,1);

log_likelihood1 = zeros(length(mu_range),1);
log_likelihood2 = zeros(length(mu_range),1);

for j = 1:obsv
    for i = 1:length(mu_range)
        mu = mu_range(i);
        log_likelihood1(i,j) = -0.5*length(data1)*log(2*pi) - length(data1)*log(sigma1)...
            - 0.5*sum(((data1(j)-mu)/sigma1).^2);
        log_likelihood2(i,j) = -0.5*length(data2)*log(2*pi) - length(data2)*log(sigma2)...
            - 0.5*sum(((data2(j)-mu)/sigma2).^2);
    end
end

f2 = figure;
t1 = tiledlayout(1,2);
t1.TileSpacing = 'compact';
t1.Padding = 'compact';

nexttile
hold on;
for j = 1:obsv
    plot(mu_range,log_likelihood1(:,j), 'Color',[0.0000 0.3176 0.6196]);
end
hold off
xlim([mu_range(1) mu_range(end)])
xlabel('Parameter $\mu$');
ylabel('log-Likelihood $\mathcal{L}$');

nexttile
hold on
for j = 1:obsv
    plot(mu_range,log_likelihood2(:,j), 'Color',[0.0000 0.3176 0.6196]);
end
hold off;
xlim([mu_range(1) mu_range(end)])
xlabel('Parameter $\mu$');
ylabel('log-Likelihood $\mathcal{L}$');


%% Figure 3.3 - Confidence Ellipses

