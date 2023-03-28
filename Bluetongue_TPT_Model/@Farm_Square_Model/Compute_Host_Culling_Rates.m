function [dairyCattleCullingRateMatrix, beefCattleCullingRateMatrix, ...
sheepCullingRateMatrix] = Compute_Host_Culling_Rates(par)

    % Computes the values of the "culling matrices" for dairy cattle, beef cattle and
    % sheep. These matrices are the same sizes as the corresponding host state matrices and
    % their (i,j)th value gives the per-capita rate of transition for individuals in the
    % (i,j)th compartment of the state matrix out of the system due to culling/slaughter.
    %
    % par: The within-farm BTV parameters struct. Bluetongue_TPT_Model.Define_Within_Farm_
    % Parameters() generates a structure of the required form.
    %
    % Returns three values: the dairy culling matrix, the beef culling matrix and the
    % sheep culling matrix.
    %
    % AUTHOR: Laurence Dhonau.

    dairyCattleCullingRateMatrix = zeros(par.L + par.M + 2, par.K + 2);

    % Culling of dairy cattle in the postpartum compartments.
    dairyCattleCullingRateMatrix(end, :) = par.chi;

    beefCattleCullingRateMatrix = zeros(par.Ltilde + par.M + 2, par.K + 2);

    % Slaughter of beef cattle in the final youth compartments.
    beefCattleCullingRateMatrix(par.Ltilde, :) = (1 - par.e) * par.Ltilde * par.omegatilde;

    % Culling of beef cattle in the postpartum compartments.
    beefCattleCullingRateMatrix(end, :) = par.chitilde;

    % Natural mortality/culling/slaughter of sheep.
    sheepCullingRateMatrix = par.b_S * ones(1, par.Kbar + 2);
end