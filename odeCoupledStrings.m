function [x1, x2] = odeCoupledStrings(t, par, con)
    
    % Check if parameters are entered
    if nargin < 3
        con = struct();        
        if nargin < 2
            par = struct();
            if nargin < 1
                error('Please insert the time variable');
            end
        end
    end
    
    % Set default values for conditions (if so)
    if ~isfield(con, 'x10'),    con.x10 = 1; end
    if ~isfield(con, 'x20'),    con.x20 = 2; end
    
    % Set default values for parameters (if so)
    if ~isfield(par, 'm1'),     par.m1  = 1; end
    if ~isfield(par, 'm2'),     par.m2  = 1; end
    if ~isfield(par, 'mk1'),    par.mk1 = 0; end
    if ~isfield(par, 'mk2'),    par.mk2 = 0; end
    if ~isfield(par, 'k1'),     par.k1  = 6; end
    if ~isfield(par, 'k2'),     par.k2  = 6; end
    if ~isfield(par, 'g'),      par.g   = 9.8; end
    
    % Total Mass [kg]
    m = par.m1 + par.m2 + par.mk1 + par.mk2;
    
    % Total Stiffness [N/m]
    k = par.k1 + par.k2;
    
    % Portion of Total Mass
    rho_m = (par.m1 + par.mk1) / m;
    
    % Portion of Total Stiffness
    rho_k = par.k1 / k;
    
    % Set Initial Conditions
    initialConditions = [con.x10, con.x20, 0, 0];
    
    % Obtain the numerical solution
    [~, sol] = ode45(@(t,X) firstOrderSystem(t,X), t, initialConditions);
    
    % Extract the variables
    x1 = sol(:, 1)';
    x2 = sol(:, 2)';
    
    % Define the model to use
    function dX = firstOrderSystem(time, X)
        dX      = nan(size(X));
        dX(1)   = X(3);
        dX(2)   = X(4);            
        dX(3)   = (-rho_k * k * X(1) + (1 - rho_k) * k * ...
            (X(2) - X(1))) / (rho_m * m);
        dX(4)   = -(1 - rho_k) * k * (X(2) - X(1)) / ((1 - rho_m) * m);       
    end
    
end