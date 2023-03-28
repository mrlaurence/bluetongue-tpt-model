function [dairyCattleMortalityRateMatrix, beefCattleMortalityRateMatrix, ...
sheepMortalityRateMatrix] = Compute_Host_Mortality_Rates(par)

    % Computes the values of the "mortality matrices" for dairy cattle, beef cattle and
    % sheep. These matrices are the same sizes as the corresponding host state matrices
    % and their (i,j)th value gives the per-capita rate of transition for individuals
    % in the (i,j)th compartment out of the system due to BTV-associated mortality.
    %
    % par: The within-farm BTV parameters struct. Bluetongue_TPT_Model.Define_Within_Farm_
    % Parameters() generates a structure of the required form.
    %
    % Returns three values: the dairy mortality matrix, the beef mortality matrix and the
    % sheep mortality matrix.
    %
    % AUTHOR: Laurence Dhonau.

    dairyCattleMortalityRateMatrix = zeros(par.L + par.M + 2, par.K + 2);

    % Disease-associated mortality of dairy cattle is assumed to occur in all infected
    % compartments at a constant rate.
    dairyCattleMortalityRateMatrix(:, 2 : end - 1) = par.d_C;

    % Disease-associated mortality of beef cattle is assumed to occur in all infected
    % compartments at a constant rate.
    beefCattleMortalityRateMatrix = zeros(par.Ltilde + par.M + 2, par.K + 2);
    beefCattleMortalityRateMatrix(:, 2 : end - 1) = par.d_C;

    % Disease-associated mortality of sheep is assumed to occur in all infected
    % compartments at a constant rate.
    sheepMortalityRateMatrix = [0, par.d_S * ones(1, par.Kbar), 0];
end