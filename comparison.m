set(0, ...
    'DefaultAxesFontSize', 20,                  ...
    'DefaultLineLineWidth', 1,                  ...
    'DefaultAxesLineWidth', 1,                  ...
    'DefaultAxesTickLabelInterpreter', 'LaTeX', ...
    'DefaultLegendInterpreter', 'LaTeX',        ...
    'DefaultFigureColor', 'White'               ...
);

%% Get experimental data

% Raw data
expData1 =  csvread('caso1.csv'); % readmatrix('caso1excel.xlsx'); %
expData2 =  csvread('caso2.csv'); % readmatrix('caso2excel.xlsx'); %

% Set initial time sample to trim signals
ts01 = 4;
ts02 = 2;

% Preprocess time data to avoid repeated values
detectValidData = @(x) ([1, diff(x(:)')] ~= 0);

validDataMask1 = detectValidData(expData1(:,1));
validDataMask1(1 : ts01) = 0;
validDataMask1(expData1(:,1) >= 10) = 0;

validDataMask2 = detectValidData(expData2(:,1));
validDataMask2(1 : ts02) = 0;
validDataMask2(expData2(:,1) >= 10) = 0;

% Define the new datasets
dataset1 = struct(          ...
    'time', expData1(validDataMask1, 1) - expData1(ts01 + 1, 1),  ...
    'x1', expData1(validDataMask1, 2),    ...
    'x2', expData1(validDataMask1, 3)     ...
    );

dataset2 = struct(          ...
    'time', expData2(validDataMask2, 1) - expData2(ts02 + 1, 1),  ...
    'x1', expData2(validDataMask2, 2),    ...
    'x2', expData2(validDataMask2, 3)     ...
    );

%% Set parameters 

% Case 1

% Set default values for conditions (if so)
con1 = struct(              ...
    't0', 0,                ...
    'x10', dataset1.x1(1),  ... % m 
    'x20', dataset1.x2(1) ... % m 
    );
    %'X1', max(dataset1.x1) - min(dataset1.x1),  ... % m 
    %'X2', max(dataset1.x2) - min(dataset1.x2) ... % m 
    %);

% Set default values for parameters (if so)
par1 = struct(          ... 
    'm1', 0.12723,      ... % kg
    'm2', 0.12781,      ... % kg
    'mk1', 0.0,      ... % kg  0.00447
    'mk2', 0.0,      ... % kg  0.00447
    'k1', 10.57372949,     ... % N/m
    'k2', 9.792093254      ... % N/m
    );

% Case 2

% Set default values for conditions (if so)
con2 = struct(              ...
    't0', 0,                ...
    'x10', dataset2.x1(1),  ... % m 
    'x20', dataset2.x2(1) ... % m 
    );
    %'X1', max(dataset2.x1) - min(dataset2.x1),  ... % m 
    %'X2', max(dataset2.x2) - min(dataset2.x2) ... % m 
    %);

% Set default values for parameters (if so)
par2 = struct(          ... 
    'm1', 0.12723,      ... % kg
    'm2', 0.20306,      ... % kg
    'mk1', 0.0,      ... % kg  0.00447
    'mk2', 0.0,      ... % kg  0.00449 
    'k1', 10.57372949,     ... % N/m
    'k2', 9.792093254      ... % N/m
    );

% Save data
save('experimental_data.mat', 'dataset1', 'dataset2', ...
    'con1', 'par1', 'con2', 'par2');

%% Plot experimental data with a sample gamma

plotModels(dataset1, par1, con1, 0.9977);

plotModels(dataset2, par2, con2, 0.9967);

%%

% t = linspace(0, 100, 5000);
% g = @(x,ep) x - ep*(x).^2;
% plot(t, sin(pi*t), t, sin(pi*g(t, 1e-4))),