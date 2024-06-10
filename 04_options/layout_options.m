function layout = layout_options

%% Appearance
colors = [0.0000 0.3176 0.6196; % dark blue
    1.0000 0.8353 0.0000; % yellow
    0.0000 0.7451 1.0000; % light blue
    0.2431 0.2667 0.2980; % dark grey
    0.6235 0.6000 0.5961]; % light grey
set(0,'defaultAxesColorOrder',colors)
set(0,'defaultAxesXGrid','on')
set(0,'defaultAxesYGrid','on')
set(0,'defaulttextinterpreter','latex');
set(0,'defaultAxesTitle','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',10)
set(0,'defaultLegendFontSize',10)
set(0,'DefaultLegendFontSizeMode','manual')
set(0, 'DefaultAxesBox', 'on');
set(0,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},...
    {'k','k','k'})

layout = struct('colors',colors);