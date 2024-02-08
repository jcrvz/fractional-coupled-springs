function [FVU, par] = objectiveFunction(dataset, par, con, model, variables)
    normaliseData = @(x) (x - min(x)) / (max(x) - min(x)) - 0.5;
    
    % Variables: [t0, x10, x20, dk1, dk2, dm1, dm2, gamma]
    if length(variables) >= 4
        if length(variables) < 7 % medium
            par.k1 = par.k1 + variables(1);
            par.k2 = par.k2 + variables(2);
            par.m1 = par.m1 + variables(3);
            par.m2 = par.m2 + variables(4);
        else % long
            par.tcoef0 = variables(1);
            par.tcoef1 = variables(2);
            par.tcoef2 = variables(3);
            par.k1 = par.k1 + variables(4);
            par.k2 = par.k2 + variables(5);
            par.m1 = par.m1 + variables(6);
            par.m2 = par.m2 + variables(7);
        end
    else % short
        par.tcoef0 = variables(1);
        par.tcoef1 = variables(2);
        par.tcoef2 = variables(3);
    end
    %     time_modified = variables(1) * dataset.time;
    
    % Choose the corresponding model and evaluate data
    if strcmp(model, 'modelCoupledSprings')
        % Get data from basic model
        [x1, x2] = modelCoupledSprings(dataset.time, par, con);
    else
        % Get data from fractional model
        [x1, x2] = fractionalCoupledSprings(dataset.time, variables(end), ...
            par, con);
        %         [x1, x2] = Double_Mass_Spring_CF2(variables(4), dataset.time, ...
        %             par.m1, par.k1, con.x10, 0, par.m2, par.k2, con.x20, 0);
    end
    
    % Compare against experimental data (least-squares)
    totalValues1 =  (normaliseData(dataset.x1) - ...
        mean(normaliseData(dataset.x1))).^2;
    totalValues2 =  (normaliseData(dataset.x2) - mean(normaliseData(dataset.x2))).^2;
    residualsValues1 = (normaliseData(dataset.x1) - normaliseData(x1)).^2;
    residualsValues2 = (normaliseData(dataset.x2) - normaliseData(x2)).^2;
    
    % Fraction of variance unexplained (FVU)
    FVU1 = sum(residualsValues1) / sum(totalValues1);
    FVU2 = sum(residualsValues2) / sum(totalValues2);
    
    % Average FVU
    FVU = (FVU1 + FVU2) / 2;
end