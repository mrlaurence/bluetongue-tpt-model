function [dairyCattleInfectionRateMatrix, beefCattleInfectionRateMatrix, ...
sheepInfectionRateMatrix] = Compute_Host_Infection_Rates(par, t, temp, H_D, H_B, H_S, ...
infectedVectorFraction)

    % Computes the values of the "infection transition matrices" for dairy cattle, beef
    % cattle and sheep. These matrices are the same sizes as the corresponding host state
    % matrices and their (i,j)th value gives the per-capita rate of transition for individuals
    % in the (i,j)th compartment of the state matrix to the compartment directly
    % rightward. Note that the rightmost columns of these transition matrices are always 0.
    %
    % par: The within-farm BTV parameters struct. Bluetongue_TPT_Model.Define_Within_Farm_
    % Parameters() generates a structure of the required form.
    % t: The current simulation time.
    % temp: The current farm square (mean daily) temperature.
    % H_D, H_B, H_S: The current total dairy cattle/beef cattle/sheep populations in the
    % farm square (not including deceased indiviudals).
    % infectedVectorFraction: The current proportion of vectors in the farm square which
    % are infected.
    %
    % Returns three values: the dairy infection transition matrix, the beef infection
    % transition matrix and the sheep infection transition matrix.
    %
    % AUTHOR: Laurence Dhonau.

    vectorActivity = Farm_Square_Model.Compute_Vector_Activity(par, t, par.activityNormaliser);

    % Compute the reciprocal of the time interval between vector blood meals.
    a = 0.0002 * temp * (temp - 3.7) * (41.9 - temp)^(1/2.7);

    % Compute the proportion of vector bites on cattle.
    phi = (H_D + H_B) / (H_D + H_B + par.sigma * H_S);

    % Compute the force of infection as a 3-element vector (giving values for dairy
    % cattle, beef cattle and sheep).
    if H_D + H_B == 0
        lambda = par.b * a * [phi, phi, (1-phi)] .* [0, 0, par.m_S] * vectorActivity ...
        * infectedVectorFraction;
    else
        lambda = par.b * a * [phi, phi, (1-phi)] .* [H_D / (H_D + H_B) * par.m_C, ...
        H_B / (H_D + H_B) * par.m_C, par.m_S] * vectorActivity * infectedVectorFraction;
    end

    dairyCattleInfectionRateMatrix = repmat([lambda(1), par.K * par.r * ...
    ones(1, par.K), 0], par.L + par.M + 2, 1);

    beefCattleInfectionRateMatrix = repmat([lambda(2), par.K * par.r * ...
    ones(1, par.K), 0], par.Ltilde + par.M + 2, 1);

    sheepInfectionRateMatrix = [lambda(3), par.Kbar * par.rbar * ones(1, par.Kbar), 0];
end