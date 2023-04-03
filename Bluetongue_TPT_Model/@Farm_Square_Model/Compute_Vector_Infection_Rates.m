function vectorInfectionRateMatrix = Compute_Vector_Infection_Rates(par, t, temp, phi, ...
infectedCattleFraction, infectedSheepFraction)

    % Computes the values of the "infection transition matrix" for vectors. This matrix is
    % the same size as the corresponding vector state matrix and its (i,j)th value gives
    % the per-capita rate of transition for individuals in the (i,j)th compartment of the
    % state matrix to the compartment directly rightward. Note that the rightmost element of
    % this transition matrix is always 0.
    %
    % par: The within-farm BTV parameters struct. Bluetongue_TPT_Model.Define_Within_Farm_
    % Parameters() generates a structure of the required form.
    % t: The current simulation time.
    % temp: The current farm square (mean daily) temperature.
    % phi: The current proportion of vector bites on cattle.
    % infectedCattleFraction, infectedSheepFraction: The current proportion of
    % cattle/sheep in the farm square which are infected.
    %
    % Returns three values: the vector infection transition matrix.
    %
    % AUTHOR: Laurence Dhonau.

    vectorActivity = Farm_Square_Model.Compute_Vector_Activity(par, t, par.activityNormaliser);

    % Compute the reciprocal of the time interval between vector blood meals.
    if temp >= 3.7 && temp <= 41.9
        a = 0.0002 * temp * (temp - 3.7) * (41.9 - temp)^(1/2.7);
    else
        a = 0;
    end

    % Compute the force of infection.
    lambda_V = par.beta * a * vectorActivity * (phi * infectedCattleFraction + ...
    (1-phi) * infectedSheepFraction);
    
    % Compute the reciprocal of the mean vector extrinsic incubation period (EIP).
    if temp >= par.T_min
        nu = par.alpha * (temp - par.T_min);
    else
        nu = 0;
    end
    
    vectorInfectionRateMatrix = [lambda_V, par.k * nu * ones(1, par.k), 0];
end