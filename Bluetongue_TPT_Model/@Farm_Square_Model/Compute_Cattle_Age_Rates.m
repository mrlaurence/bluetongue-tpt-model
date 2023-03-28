function [dairyCattleAgeRateMatrix, beefCattleAgeRateMatrix] = Compute_Cattle_Age_Rates(par)

    % Computes the values of the "age transition matrices" for dairy and beef cattle.
    % These matrices are the same size as the corresponding dairy/beef state matrix and
    % their (i,j)th value gives the per-capita rate of transition for individuals in the
    % (i,j)th compartment of the state matrix to the compartment directly below (or, in
    % the case of the postpartum compartments, back to the corresponding adult compartment).
    %
    % par: The within-farm BTV parameters struct. Bluetongue_TPT_Model.Define_Within_Farm_
    % Parameters() generates a structure of the required form.
    %
    % Returns two values: the dairy age transition matrix and the beef age transition
    % matrix.
    %
    % AUTHOR: Laurence Dhonau.

    dairyCattleAgeRateMatrix = zeros(par.L + par.M + 2, par.K + 2);

    % Youth compartment -> next youth compartment (or adult compartment).
    dairyCattleAgeRateMatrix(1 : par.L, :) = par.L * par.omega;

    % Adult compartment -> 1st pregnancy compartment.
    dairyCattleAgeRateMatrix(par.L + 1, :) = par.psi;

    % Pregnancy compartment -> next pregnancy compartment (or postpartum compartment).
    dairyCattleAgeRateMatrix(par.L + 2 : par.L + par.M + 1, :) = par.M * par.delta;

    % Postpartum compartment -> adult compartment.
    dairyCattleAgeRateMatrix(end, :) = par.epsilon - par.chi;

    beefCattleAgeRateMatrix = zeros(par.Ltilde + par.M + 2, par.K + 2);

    % Youth compartment -> next youth compartment. 
    beefCattleAgeRateMatrix(1 : par.Ltilde - 1, :) = par.Ltilde * par.omegatilde;

    % Final youth compartment -> adult compartment.
    beefCattleAgeRateMatrix(par.Ltilde, :) = par.e * par.Ltilde * par.omegatilde;

    % Adult compartment -> 1st pregnancy compartment.
    beefCattleAgeRateMatrix(par.Ltilde + 1, :) = par.psi;

    % Pregnancy compartment -> next pregnancy compartment (or postpartum compartment).
    beefCattleAgeRateMatrix(par.Ltilde + 2 : par.Ltilde + par.M + 1, :) = par.M * par.delta;

    % Postpartum compartment -> adult compartment.
    beefCattleAgeRateMatrix(end, :) = par.epsilontilde - par.chitilde;
end