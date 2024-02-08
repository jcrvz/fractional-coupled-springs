function [x1, x2] = modelCoupledSprings(t, par, con)
    
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
    if ~isfield(con, 'x10'),    con.x10 = nan; end
    if ~isfield(con, 'x20'),    con.x20 = nan; end
    if ~isfield(con, 't0'),     con.t0 = 0; end
    
    % Set default values for parameters (if so)
    if ~isfield(par, 'm1'),     par.m1  = nan; end
    if ~isfield(par, 'm2'),     par.m2  = nan; end
    if ~isfield(par, 'mk1'),    par.mk1 = nan; end
    if ~isfield(par, 'mk2'),    par.mk2 = nan; end
    if ~isfield(par, 'k1'),     par.k1  = nan; end
    if ~isfield(par, 'k2'),     par.k2  = nan; end
    if ~isfield(par, 'tcoef0'), par.tcoef0  = 0; end
    if ~isfield(par, 'tcoef1'), par.tcoef1  = 1; end
    if ~isfield(par, 'tcoef2'), par.tcoef2  = 0; end
    
    % Total Mass [kg]
    m = par.m1 + par.m2 + par.mk1 + par.mk2;
    
    % Total Stiffness [N/m]
    k = par.k1 + par.k2;
    
    % Portion of Total Mass
    rho_m = (par.m1 + par.mk1) / m;
    
    % Portion of Total Stiffness
    rho_k = par.k1 / k;
    
    % Dimensional length scale [m]
    X = 1 / k;
    if ~isfield(con, 'X1'),    con.X1 = X; end
    if ~isfield(con, 'X2'),    con.X2 = X; end
    
    % Dimensional temporal scale [s]
    T = sqrt(m / k);
    
    % Determine the non-dimensional initial conditions
    x10_nd = con.x10 / con.X1;
    x20_nd = con.x20 / con.X2;
    
    % Find the frequencies
    alpha_part_1 = (1 - rho_k * rho_m) / (2 * rho_m * (1 - rho_m));
    alpha_part_2 = 4 * rho_k * rho_m * (1 - rho_k) * (1 - rho_m) / ...
        (1 - rho_k * rho_m)^2;
    alpha2_1 = alpha_part_1 * (1 - sqrt(1 - alpha_part_2));
    alpha2_2 = alpha_part_1 * (1 + sqrt(1 - alpha_part_2));
    
%     alphas = roots([1, (1 - rho_k * rho_m) / (rho_m * (1 - rho_m)), ...
%     rho_k * (1 - rho_k) / (rho_m * (1 - rho_m))]);
%     alpha2_1 = -alphas(1);
%     alpha2_2 = -alphas(2);
    
    % Coefficients
    A0 = 1 / ((alpha2_1 - alpha2_2) * (1 - rho_m) * rho_m);
    A1 = alpha2_1 * rho_m * (1 - rho_m) * x10_nd - ...
        (1 - rho_k) * (rho_m * x10_nd + (1 - rho_m) * x20_nd);
    A2 = - alpha2_2 * rho_m * (1 - rho_m) * x10_nd + ...
        (1 - rho_k) * (rho_m * x10_nd + (1 - rho_m) * x20_nd);
    A3 = alpha2_1  * rho_m * (1 - rho_m) * x20_nd - ...
        (1 - rho_m) * x20_nd - rho_m * (1 - rho_k) * x10_nd;
    A4 = - alpha2_2 * rho_m * (1 - rho_m) * x20_nd + ...
        (1 - rho_m) * x20_nd + rho_m * (1 - rho_k) * x10_nd;
    
    % Get the non-dimensional variable t
    t_nd = (par.tcoef0 + par.tcoef1 * t + par.tcoef2 * t.^2) / T;
    
    % Find the non-dimensional variables x1 and x2    
    x1_nd = A0 * (  A1 * cos(sqrt(alpha2_1) * t_nd) + ...
                    A2 * cos(sqrt(alpha2_2) * t_nd));
%         - (1 - rho_k) * sin(sqrt(alpha2_1) * t_nd) / sqrt(alpha2_1) ...
%         + (1 - rho_k) * sin(sqrt(alpha2_2) * t_nd) / sqrt(alpha2_2) ...
    x2_nd = A0 * (  A3 * cos(sqrt(alpha2_1) * t_nd) + ...
                    A4 * cos(sqrt(alpha2_2) * t_nd));
%         1 / (alpha2_1 * alpha2_2);
    
    % Find the dimensional variables x1 and x2
    x1 = x1_nd * con.X1 ;
    x2 = x2_nd * con.X2 ;    
end