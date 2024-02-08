function [x1, x2] = fractionalCoupledSprings(t, gamma, par, con)
    
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
    
    % Parameter for determinants
    u_1 = (1 - rho_k) *  (rho_m * x10_nd + (1 - rho_m) * x20_nd);
    u_2 = rho_m * (1 - rho_k) * x10_nd + (1 - rho_m) * x20_nd;
    
    % Find the fractional damping factor and frequency
    beta_1 = (1 - gamma) * alpha2_1 / 2;
    beta_2 = (1 - gamma) * alpha2_2 / 2;
    omega_1 = sqrt(gamma * alpha2_1 - (1 - gamma)^2 * alpha2_1^2 / 4);
    omega_2 = sqrt(gamma * alpha2_2 - (1 - gamma)^2 * alpha2_2^2 / 4);
    
    % Coefficients
    coefficient_matrix = [
        1, 0, 1, 0;
        2*beta_2, 1, 2*beta_1, 1;
        beta_2^2+omega_2^2, 2*beta_2, beta_1^2+omega_1^2, 2*beta_1;
        0, beta_2^2+omega_2^2, 0, beta_1^2+omega_1^2
        ];
    coeff1 = coefficient_matrix \ [ ...
        rho_m*(1-rho_m)*x10_nd; (1-gamma)*u_1; gamma*u_1; 0];
    coeff2 = coefficient_matrix \ [ ...
        rho_m*(1-rho_m)*x20_nd; (1-gamma)*u_2; gamma*u_2; 0];
    
    % Overall scale
    A0 = 1 / (rho_m * (1 - rho_m));
    
    % Get the non-dimensional variable t
    t_nd = (par.tcoef0 + par.tcoef1 * t + par.tcoef2 * t.^2) / T;
    
    % Find the non-dimensional variables x1 and x2
    x1_nd = A0 * ( exp(-beta_1 * t_nd) .* ( ...
        coeff1(1) * cos(omega_1 * t_nd) + ...
        (coeff1(2) - beta_1 * coeff1(1)) * ...
        sin(omega_1 * t_nd) / omega_1 ) + ...
        exp(-beta_2 * t_nd) .* ( ...
        coeff1(3) * cos(omega_2 * t_nd) + ...
        (coeff1(4) - beta_2 * coeff1(3)) * ...
        sin(omega_2 * t_nd) / omega_2));
    x2_nd = A0 * ( exp(-beta_1 * t_nd) .* ( ...
        coeff2(1) * cos(omega_1 * t_nd) + ...
        (coeff2(2) - beta_1 * coeff2(1)) * ...
        sin(omega_1 * t_nd) / omega_1 ) + ...
        exp(-beta_2 * t_nd) .* ( ...
        coeff2(3) * cos(omega_2 * t_nd) + ...
        (coeff2(4) - beta_2 * coeff2(3)) * ...
        sin(omega_2 * t_nd) / omega_2));
    
    % Find the dimensioal variables x1 and x2
    x1 = x1_nd * con.X1;
    x2 = x2_nd * con.X2;
end