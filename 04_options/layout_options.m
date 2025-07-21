function layout = layout_options

%% Appearance
colors = [0.0000 0.3176 0.6196; % dark blue
    0.9569 0.5961 0.6118; % salmon pink
    0.4000 0.8000 0.2000 % green
    0.9216 0.8235 0.705; % bisque
    0.6235 0.6000 0.5961]; % light grey
set(0,'defaultAxesColorOrder',colors)
set(0,'defaultAxesXGrid','on')
set(0,'defaultAxesYGrid','on')
set(0,'defaulttextinterpreter','latex');
set(0,'defaultAxesTitle','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',9)
set(0,'defaultLegendFontSize',9)
set(0,'DefaultLegendFontSizeMode','manual')
set(0, 'DefaultAxesBox', 'on');
set(0,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},...
    {'k','k','k'})
set(0,'DefaultFigureColormap',feval('sky'));

layout = struct('colors',colors);