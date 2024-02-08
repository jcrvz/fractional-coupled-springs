set(0, ...
    'DefaultAxesFontSize', 20,                  ...
    'DefaultLineLineWidth', 1,                  ...
    'DefaultAxesLineWidth', 1,                  ...
    'DefaultAxesTickLabelInterpreter', 'LaTeX', ...
    'DefaultLegendInterpreter', 'LaTeX',        ...
    'DefaultFigureColor', 'White'               ...
    );

% Time parameters
time = linspace(0, 10, 1000);

% Set default values for conditions (if so)
con = struct(           ...
    'x10', 0.01,        ... % m
    'x20', 0.01         ... % m
    );

% Set default values for parameters (if so)
par = struct(           ...
    'm1', 0.100,        ... % kg
    'm2', 0.100,        ... % kg
    'mk1', 0,           ... % kg  0.00447
    'mk2', 0,           ... % kg  0.00449
    'k1', 10,     ... % N/m
    'k2', 10      ... % N/m
    );

%%

% Get Simulink model
simOut = sim('simCoupledStrings_R2016a.slx', time);
data = get(simOut, 'data');

% Get data from basic model
[x1ord, x2ord] = modelCoupledSprings(time, par, con);

% Set the fractional order
gamma = 0.99;

% Get data from basic model
[x1f, x2f] = fractionalCoupledSprings(time, gamma, par, con);

%% Comparative between models
% saveas(get_param('simCoupledStrings','Handle'), 'model', 'svg')

% Prepare the figures results
Fi1 = figure('Name', 'Models Comparison 1', ...
    'Units', 'normalized', 'Position', [0.2677 0.4574 0.4703 0.2657]);
Ax1 = axes('NextPlot', 'Add', 'Box', 'On');

plot(Ax1, data.Time, data.Data(:,1), '-r', 'LineWidth', 2, ...
    'DisplayName', 'Simulink');
plot(Ax1, time, x1ord, '--g', ...
    'DisplayName', 'Ordinary');
plot(Ax1, time, x1f, '-.b', ...
    'DisplayName', sprintf('CF $$\\gamma = %.2f$$', gamma));

xlabel('Time, $$t$$ [s]', 'Interpreter', 'LaTeX');
ylabel('Position $$x_1$$ [m]', 'Interpreter', 'LaTeX');
ylim([-0.02 0.02])

% Set legends
Le1 = legend(Ax1, 'Show', 'Orientation', 'Horizontal', 'Box', 'Off', ...
    'Location', 'North');

%% Ordinary vs. Fractional

% Prepare the figures results
Fi2 = figure('Name', 'Models Comparison 2', ...
    'Units', 'normalized', 'Position', [0.2677 0.1231 0.4703 0.2657]);
Ax2 = axes('NextPlot', 'Add', 'Box', 'On');

plot(Ax2, data.Time, data.Data(:,2), '-r', 'LineWidth', 2, ...
    'DisplayName', 'Simulink');
plot(Ax2, time, x2ord, '--g', ...
    'DisplayName', 'Ordinary');
plot(Ax2, time, x2f, '-.b', ...
    'DisplayName', sprintf('CF $$\\gamma = %.2f$$', gamma));

xlabel('Time, $$t$$ [s]', 'Interpreter', 'LaTeX');
ylabel('Position $$x_2$$ [m]', 'Interpreter', 'LaTeX');
ylim([-0.02 0.02])

% Set legends
%Le2 = legend(Ax2, 'Show', 'Orientation', 'Horizontal');

%% Variation of gamma

n_gamma = 4;

gamma_values = linspace(0.95, 1.0, n_gamma);

X1f = nan(n_gamma, 1000); X2f = nan(n_gamma, 1000);

[G, T] = meshgrid(time, gamma_values);

% Prepare the figures results
Fi3 = figure('Name', 'Models Comparison 3', ...
    'Units', 'normalized', 'Position', [0.1153 0.2256 0.3521 0.5300]);
Ax3 = axes('NextPlot', 'Add', 'Box', 'On');

% Prepare the figures results
Fi4 = figure('Name', 'Models Comparison 4', ...
    'Units', 'normalized', 'Position', [0.5188 0.2078 0.3521 0.5300]);
Ax4 = axes('NextPlot', 'Add', 'Box', 'On');

colormap jet

for ii = 1 : n_gamma
    
    gamma = gamma_values(ii);
    
    % Get data from basic model
    [x1f, x2f] = fractionalCoupledSprings(time, gamma, par, con);
    
    X1f(ii, :) = x1f;
    X2f(ii, :) = x2f;
    
    plot3(Ax3, time, gamma * ones(1, 1000), x1f);
    plot3(Ax4, time, gamma * ones(1, 1000), x2f);
    
end

view(Ax3, 46, 36);  view(Ax4, 46, 36);
set([Ax3, Ax4], 'YTick', gamma_values, 'YGrid', 'on', ...
    'ZLim', [-0.02 0.02]);
xlabel([Ax3, Ax4], 'Time, $$t$$ [s]', 'Interpreter', 'LaTeX');
ylabel([Ax3, Ax4], 'Fractional order $$\gamma$$','Interpreter', 'LaTeX');
zlabel(Ax3, 'Position $$x_1$$ [m]', 'Interpreter', 'LaTeX');
zlabel(Ax4, 'Position $$x_2$$ [m]', 'Interpreter', 'LaTeX');
