set(0, ...
    'DefaultAxesFontSize', 20,                  ...
    'DefaultLineLineWidth', 1,                  ...
    'DefaultAxesLineWidth', 1,                  ...
    'DefaultAxesTickLabelInterpreter', 'LaTeX', ...
    'DefaultLegendInterpreter', 'LaTeX',        ...
    'DefaultFigureColor', 'White'               ...
    );


%% Get experimental data

load('experimental_data.mat')

runOptimisation = true;
% Set case
expCase     = 1;    % 1 or 2
modeParams  = '_3m_m';  % s or l

% Get the parameters
dataset = eval(sprintf('dataset%d', expCase));
par = eval(sprintf('par%d', expCase));
con = eval(sprintf('con%d', expCase));


fileName = sprintf('FRcase%d%s', expCase, modeParams);

%% Set the objective function

% Available models
models = {'modelCoupledSprings', 'fractionalGarrapaCoupledSprings', ...
    'fractionalCoupledSprings'};
prettyNames = {'Ordinary', 'Caputo', 'Caputo-Fabrizio'};

if runOptimisation
    % Get problems
    problem1 = @(vars) objectiveFunction(dataset, par, con, models{1}, vars);
    problem2 = @(vars) objectiveFunction(dataset, par, con, models{2}, vars);
    problem3 = @(vars) objectiveFunction(dataset, par, con, models{3}, vars);
    
    if fileName(end) == 's' % Short
        lb = [0 0.5 -0.1];
        ub = [0.5 1.0 0.1];
    elseif fileName(end) == 'm' % Medium
        % variables = [dt0, dx10, dx20, dk1, dk2, dm1, dm2, gamma]
        p = 0.5;
        lb = [0, 0, 0, 0];
        ub = [p*par.k1, p*par.k2, p*par.m1, p*par.m2];
    else % Long
        % variables = [dt0, dx10, dx20, dk1, dk2, dm1, dm2, gamma]
        p = 0.5;
        lb = [0 0.5 -0.1 0, 0, 0, 0];
        ub = [0.5 1.0 0.1 p*par.k1, p*par.k2, p*par.m1, p*par.m2];
    end
    numOfRepetitions = 50;
    
    DataXSol1 = nan(numOfRepetitions, length(lb));
    DataFX1 = nan(numOfRepetitions, 1);
    
    DataXSol2 = nan(numOfRepetitions, length(lb)+1);
    DataFX2 = nan(numOfRepetitions, 1);
    
    DataXSol3 = nan(numOfRepetitions, length(lb)+1);
    DataFX3 = nan(numOfRepetitions, 1);
    
    for ii = 1 : numOfRepetitions
        [xSol1_, fx1_] = CSOA(problem1, [lb', ub']);
        [xSol2_, fx2_] = CSOA(problem2, [[lb, 0.9]',[ub, 1.0-eps]']);
        [xSol3_, fx3_] = CSOA(problem3, [[lb, 0.9]',[ub, 1.0-eps]']);
        
        mark1 = ' '; mark2 = ' '; mark3 = ' ';
        if ii == 1
            xSol1 = xSol1_; fx1 = fx1_;
            xSol2 = xSol2_; fx2 = fx2_;
            xSol3 = xSol3_; fx3 = fx3_;
        else
            if fx1_ < fx1,  xSol1 = xSol1_; fx1 = fx1_;  mark1 = '*'; end
            if fx2_ < fx2,  xSol2 = xSol2_; fx2 = fx2_;  mark2 = '*'; end
            if fx3_ < fx3,  xSol3 = xSol3_; fx3 = fx3_;  mark3 = '*'; end
        end
        
        fprintf('[%d]-%s fx1: %2g, x1: [%s]\n', ii, mark1, fx1_, ...
            sprintf('%.4g ', xSol1_));
        fprintf('[%d]-%s fx2: %2g, x2: [%s]\n', ii, mark2, fx2_, ...
            sprintf('%.4g ', xSol2_));
        fprintf('[%d]-%s fx3: %2g, x3: [%s]\n', ii, mark3, fx3_, ...
            sprintf('%.4g ', xSol3_));
        
        % Save data
        DataXSol1(ii,:) = xSol1_; DataFX1(ii) = fx1_;
        DataXSol2(ii,:) = xSol2_; DataFX2(ii) = fx2_;
        DataXSol2(ii,:) = xSol3_; DataFX3(ii) = fx3_;
        
    end
    
    % save data
    save([fileName, '.mat'], ...
        'DataXSol1', 'DataXSol2', 'DataFX1', 'DataFX2', ...
        'DataFX3', 'DataFX3');
    
else
    % load data
    load([fileName, '.mat']);
    
    % Find the best solution
    [fx1, argFx1] = min(DataFX1);
    [fx2, argFx2] = min(DataFX2);
    [fx3, argFx3] = min(DataFX3);
    
    xSol1 = DataXSol1(argFx1, :);
    xSol2 = DataXSol2(argFx2, :);
    xSol3 = DataXSol2(argFx3, :);
end

%% Violin Plots
Fi0 = figure('Name', ['Violinplot', fileName]);
set(Fi0, 'Units', 'normalized', 'Position', [0.3953 0.3611 0.2089 0.3472]);
Ax0 = axes('NextPlot', 'Add', 'Box', 'On');

Vp0 = violinplot([DataFX1, DataFX2, DataFX3]);

colours = repmat(lines(4), 3, 1);
coloursCell = num2cell(colours, 2);
[Vp0.ViolinColor] = coloursCell{:};

xlabel('Model', 'Interpreter', 'LaTeX');
ylabel('FVU', 'Interpreter', 'LaTeX');
Ax0.XTickLabel = prettyNames;

%% Evaluate and plot the obtained results
normaliseData = @(x) (x - min(x)) / (max(x) - min(x)) - 0.5;

% Prepare the figures results
Fi1 = figure('Name', 'Values for x1');
set(Fi1, 'Units', 'normalized', ...
    'Position', [0.2677 0.4574 0.4703 0.2657]);
Ax1 = axes('NextPlot', 'Add', 'Box', 'On');

Fi2 = figure('Name', 'Values for x2');
set(Fi2, 'Units', 'normalized', ...
    'Position', [0.2677 0.1231 0.4703 0.2657]);
Ax2 = axes('NextPlot', 'Add', 'Box', 'On');

Fi3 = figure('Name', 'Cartesian graph');
set(Fi3, 'Units', 'normalized', ...
    'Position', [0.3948 0.1796 0.2604 0.4537]);
Ax3 = axes('NextPlot', 'Add', 'Box', 'On');

% Plot experimental data
plot(Ax1, dataset.time, normaliseData(dataset.x1), '--k', ...
    'DisplayName', 'Data');
plot(Ax2, dataset.time, normaliseData(dataset.x2), '--k', ...
    'DisplayName', 'Data');
plot(Ax3, normaliseData(dataset.x1), normaliseData(dataset.x2), ...
    '-k', 'DisplayName', 'Data')

for modelId = 1 : numel(models)
    con_mod = con;
    par_mod = par;
    
    % Set the 'optimal' variables
    %     con_mod.t0 = con_mod.t0 + eval(sprintf('xSol%d(1)', modelId));
    %     con_mod.x10 = con_mod.x10 + eval(sprintf('xSol%d(1)', modelId));
    %     con_mod.x20 = con_mod.x20 + eval(sprintf('xSol%d(2)', modelId));
    if fileName(end) == 's'
        par_mod.tcoef0 = eval(sprintf('xSol%d(1)', modelId));
        par_mod.tcoef1 = eval(sprintf('xSol%d(2)', modelId));
        par_mod.tcoef2 = eval(sprintf('xSol%d(3)', modelId));
        
    elseif fileName(end) == 'm'
        par_mod.k1 = par_mod.k1 + eval(sprintf('xSol%d(1)', modelId));
        par_mod.k2 = par_mod.k2 + eval(sprintf('xSol%d(2)', modelId));
        par_mod.m1 = par_mod.m1 + eval(sprintf('xSol%d(3)', modelId));
        par_mod.m2 = par_mod.m2 + eval(sprintf('xSol%d(4)', modelId));
        
        eval(sprintf('xSol%d(1) = par_mod.k1;', modelId));
        eval(sprintf('xSol%d(2) = par_mod.k2;', modelId));
        eval(sprintf('xSol%d(3) = par_mod.m1;', modelId));
        eval(sprintf('xSol%d(4) = par_mod.m2;', modelId));
        
    elseif fileName(end) == 'l'
        par_mod.tcoef0 = eval(sprintf('xSol%d(1)', modelId));
        par_mod.tcoef1 = eval(sprintf('xSol%d(2)', modelId));
        par_mod.tcoef2 = eval(sprintf('xSol%d(3)', modelId));
        par_mod.k1 = par_mod.k1 + eval(sprintf('xSol%d(4)', modelId));
        par_mod.k2 = par_mod.k2 + eval(sprintf('xSol%d(5)', modelId));
        par_mod.m1 = par_mod.m1 + eval(sprintf('xSol%d(6)', modelId));
        par_mod.m2 = par_mod.m2 + eval(sprintf('xSol%d(7)', modelId));
        
        eval(sprintf('xSol%d(4) = par_mod.k1;', modelId));
        eval(sprintf('xSol%d(5) = par_mod.k2;', modelId));
        eval(sprintf('xSol%d(6) = par_mod.m1;', modelId));
        eval(sprintf('xSol%d(7) = par_mod.m2;', modelId));
    end
    
    % Evaluate the model
    switch modelId
        case 1
            [x1, x2] = modelCoupledSprings(dataset.time, par_mod, con_mod);
            
            sSol1 = sprintf(' %4g &', xSol1);
            fprintf('%.4g &%s \\\\ \n', fx1, sSol1(1:end-1));
            
            lineColor = 'b';
        case 2
            gamma = eval(sprintf('xSol%d(end)', modelId));
            [x1, x2] = fractionalCoupledSprings(...
                dataset.time, gamma, par_mod, con_mod);
            
            % Print results
            sSol2 = sprintf(' %4g &', xSol2);
            fprintf('%.4g &%s \\\\ \n', fx2, sSol2(1:end-1));
            
            lineColor = 'r';
        case 3
            gamma = eval(sprintf('xSol%d(end)', modelId));
            [x1, x2] = fractionalGarrapaCoupledSprings(...
                dataset.time, gamma, par_mod, con_mod);
            
            % Print results
            sSol3 = sprintf(' %4g &', xSol3);
            fprintf('%.4g &%s \\\\ \n', fx3, sSol3(1:end-1));
            
            lineColor = 'r';
    end
    
    % Plot data from the current model
    plot(Ax1, dataset.time, normaliseData(x1), ['-', lineColor], ...
        'DisplayName', prettyNames{modelId}, 'LineWidth', 1.0);
    plot(Ax2, dataset.time, normaliseData(x2), ['-', lineColor], ...
        'DisplayName', prettyNames{modelId}, 'LineWidth', 1.0);
    plot(Ax3, normaliseData(x1), normaliseData(x2), ['-', lineColor], ...
        'DisplayName', prettyNames{modelId}, 'LineWidth', 1.0);
end

% Final adjustments
legend(Ax1, 'Show', 'Orientation', 'horizontal');
legend(Ax3, 'Show', 'Location', 'Best', 'Box', 'off');

xlabel([Ax1, Ax2], 'Time $$t$$ [s]', 'Interpreter', 'LaTeX');
xlabel(Ax3, 'Norm. Position $$\tilde{x}_1$$', 'Interpreter', 'LaTeX');
ylabel(Ax1, 'Norm. Position $$\tilde{x}_1$$', 'Interpreter', 'LaTeX');
ylabel([Ax2, Ax3], 'Norm. Position $$\tilde{x}_2$$', 'Interpreter', 'LaTeX');
xlim([Ax1, Ax2], [0, 10]); %dataset.time(end)]);
axis(Ax3, 'square')