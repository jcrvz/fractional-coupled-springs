function [x1, x2] = fractionalGarrapaCoupledSprings(t, gamma, par, con)
    
    % Check if parameters are entered
    if nargin < 4
        con = struct();
        if nargin < 3
            par = struct();
            if nargin < 2
                error('Please insert the time variable');
            end
        end
    end
    
    % Set default values for conditions (if so)
    if ~isfield(con, 'x10'),    con.x10 = nan; end
    if ~isfield(con, 'x20'),    con.x20 = nan; end
    if ~isfield(con, 't0'),     con.t0 = 0;    end
    
    % Set default values for parameters (if so)
    if ~isfield(par, 'm1'),     par.m1  = nan; end
    if ~isfield(par, 'm2'),     par.m2  = nan; end
    if ~isfield(par, 'mk1'),    par.mk1 = nan; end
    if ~isfield(par, 'mk2'),    par.mk2 = nan; end
    if ~isfield(par, 'k1'),     par.k1  = nan; end
    if ~isfield(par, 'k2'),     par.k2  = nan; end
    
    gammaVector = gamma * ones(1, 4);
    
    % Planck's time
    tp = 5.39106e-44;
    
    A = -(tp^(1 - gamma)) * (par.k1 + par.k2) / par.m1; 
    B = (tp^(1 - gamma)) * par.k2 / par.m1;
    C = (tp^(1 - gamma)) * par.k2 / par.m2; 
    D = -(tp^(1 - gamma))* par.k2 / par.m2;
    
    parameterMatrix = [A, B, C, D];
    lambdaVector    = [0 , 1, 0, 1];
    
    
    f_fun=@(r,y,par)[y(2); A*y(1) + B*y(3); y(4); C*y(1) + D*y(3)];    
    %J_fun=@(r,y,par)[0,1,0,0; A,0,B,0; 0,0,0,1; C,0,D,0];
    
    % Define~ parameters for the numerical method
    [~, y] = MT_FDE_PI1_Ex(...
        gammaVector, lambdaVector, f_fun, con.t0, t(end), ...
        [con.x10; 0; con.x20; 0], 0.0001, parameterMatrix);
    
    % Get results
    x1 = y(1,:);
    %zv = y(2,:);
    x2 = y(3,:);
    %wv = y(4,:);
end