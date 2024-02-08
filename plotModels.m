function handles = plotModels(dataset, par, con, gamma)
    normaliseData = @(x) x; %(x - min(x)) / (max(x) - min(x)) - 0.5;
    
    % Get data from basic model
    [x1ord, x2ord] = modelCoupledSprings(dataset.time, par, con);
    
    % Get data from fractional model
    [x1cfd, x2cfd] = fractionalCoupledSprings(dataset.time, gamma, par, con);
%     [x1cfd, x2cfd] = Double_Mass_Spring_CF2(gamma, dataset.time, ...
%         par.m1, par.k1, con.x10, 0, par.m2, par.k2, con.x20, 0);
    
    % ------------------------------------------------------------------------
    
    % Prepare the figures results
    Fi1 = figure('Name', 'Values for x1');
    set(Fi1, 'Units', 'normalized', ...
        'Position', [0.2677 0.4574 0.4703 0.2657]);
    Ax1 = axes('NextPlot', 'Add', 'Box', 'On');
    
    % Plot experimental data for x1
    plot(dataset.time, normaliseData(dataset.x1), 'k', ...
        'DisplayName', 'Data', 'LineWidth', 1.1);
    
    % Plot ordinary model for x1
    plot(dataset.time, normaliseData(x1ord), 'b-.', ...
        'DisplayName', 'Ordinary', 'LineWidth', 1.1);
    
    % Plot fractional model for x1
    plot(dataset.time, normaliseData(x1cfd), 'r--', 'DisplayName', ...
        sprintf('CF $$\\gamma=%.4f$$', gamma), 'LineWidth', 1.1);
    
    % Legend
    Le1 = legend('Show');
    
    % Labels
    xlabel('Time, $$t$$ [s]', 'Interpreter', 'LaTeX');
    ylabel('Position $$m_1$$, $$x_1$$ [m]', 'Interpreter', 'LaTeX');
%     xlim([0, 20])
    
    % ------------------------------------------------------------------------
    
    % Prepare the figures results
    Fi2 = figure('Name', 'Values for x2');
    set(Fi2, 'Units', 'normalized', ...
        'Position', [0.2677 0.1231 0.4703 0.2657]);
    Ax2 = axes('NextPlot', 'Add', 'Box', 'On');
    
    % Plot experimental data for x2
    plot(dataset.time, normaliseData(dataset.x2), 'k', ...
        'DisplayName', 'Data', 'LineWidth', 1.1)
    
    % Plot ordinary model for x2
    plot(dataset.time, normaliseData(x2ord), 'b-.', ...
        'DisplayName', 'Ordinary', 'LineWidth', 1.1);
    
    % Plot fractional model for x2
    plot(dataset.time, normaliseData(x2cfd), 'r--', 'DisplayName', ...
        sprintf('CF $$\\gamma=%.4f$$', gamma), 'LineWidth', 1.1);
    
    % Legend
    Le2 = legend('Show');
    
    % Labels
    xlabel('Time, $$t$$ [s]', 'Interpreter', 'LaTeX');
    ylabel('Position $$m_2$$, $$x_2$$ [m]', 'Interpreter', 'LaTeX');
%     xlim([0, 20])
    
    
    % ---
    handles = [Fi1, Ax1, Le1, Fi2, Ax2, Le2];
    
end