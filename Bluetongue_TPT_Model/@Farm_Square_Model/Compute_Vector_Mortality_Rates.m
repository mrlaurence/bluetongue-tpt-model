function vectorMortalityRateMatrix = Compute_Vector_Mortality_Rates(par, temp)

    % Computes the value of the "mortality matrix" for vectors. This matrix is the same
    % size as the vector state matrix and its (i,j)th value gives the per-capita rate of
    % transition for individuals in the (i,j)th compartment out of the system due to
    % natural mortality.
    %
    % par: The within-farm BTV parameters struct. Bluetongue_TPT_Model.Define_Within_Farm_
    % Parameters() generates a structure of the required form.
    % temp: The current farm square (mean daily) temperature.
    %
    % Returns the vector mortality matrix.

    % Compute the vector mortality rate.
    %
    % AUTHOR: Laurence Dhonau.
    
    rho = 0.009 * exp(0.16 * temp);

    vectorMortalityRateMatrix = [0, rho * ones(1, par.k + 1)];
end