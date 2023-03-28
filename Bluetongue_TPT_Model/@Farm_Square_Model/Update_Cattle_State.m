function cattleStateMatrix = Update_Cattle_State(self, L, M, cattleStateMatrix, infectionStateChange, ...
ageStateChange, mortalityStateChange, cullingStateChange, ageReEntryRow, femaleBirthsOnly)

    % Updates a dairy/beef cattle state matrix (i.e. one element of self.dairyCattleStates/
    % self.beefCattleStates) based on the passed infectionStateChange, ageStateChange,
    % mortalityStateChange and cullingStateChange matrices. For details of the effects of
    % these matrices on the cattle state, see Farm_Square_Model.Compute_Cattle_Age_Rates(),
    % Farm_Square_Model.Compute_Host_Infection_Rates(), Farm_Square_Model.Compute_Host_
    % Mortality_Rates() and Farm_Square_Model.Compute_Host_Culling_Rates(). This function
    % also incorporates the birth events described in the project report, with the
    % possibility of transplacental transmission.
    %
    % L: The number of youth stages in dairy/beef cattle.
    % M: The number of pregnancy stages in dairy/beef cattle.
    %
    % cattleStateMatrix: A matrix which spans one index (in the third dimension) of a
    % Labelled_Tensor object (either self.dairyCattleStates or self.beefCattleStates).
    %
    % infectionStateChange: A matrix comprising the number of individuals in each
    % compartment in the state matrix who should transition rightwards.
    % ageStateChange: A matrix comprising the number of individuals in each compartment in
    % the state matrix who should transition downwards.
    % mortalityStateChange: A matrix comprising the number of individuals in each
    % compartment in the state matrix who experience disease-associated mortality.
    % cullingStateChange: A matrix comprising the number of indiviudals in each
    % compartment in the state matrix who experience culling.
    % 
    % ageReEntryRow: The index of the "adult" row in the state matrix.
    % femaleBirthsOnly: A bool indicating whether only female births should enter the
    % population structure.
    %
    % Returns an updated cattle state matrix.
    %
    % AUTHOR: Laurence Dhonau.

    % Retrieve the within-farm parameters structure.
    par = self.parameters;

    % Update the cattle state matrix to incorporate infection events in this time step.
    cattleStateMatrix = cattleStateMatrix + [zeros(L + M + 2, 1), infectionStateChange(:, ...
    1:end-1)] - infectionStateChange;
    
    % Update the cattle state matrix to incorporate aging events in this time step
    % (excluding postpartum -> adult events).
    cattleStateMatrix = cattleStateMatrix + [zeros(1, par.K + 2); ageStateChange(1:end-1, ...
    :)] - ageStateChange;

    % Update the cattle state matrix to incorporate "postpartum -> adult" aging events in
    % this time step.
    cattleStateMatrix(ageReEntryRow, :) = cattleStateMatrix(ageReEntryRow, :) + ...
    ageStateChange(end,:);
    
    % Update the cattle state matrix to incorporate disease-associated mortality and
    % culling events in this time step.
    cattleStateMatrix = cattleStateMatrix - mortalityStateChange - cullingStateChange;

    %% ----- SECTION: BIRTH EVENTS -----
    
    % Add up the number of susceptible/recovered cows who exited the final pregnancy stage
    % this time step.
    birthEventsFromNonInfectedCows = ageStateChange(end - 1, 1) + ageStateChange(end - 1, end);

    % Add up the number of infected cows who exited the final pregnancy stage this time
    % step.
    birthEventsFromInfectedCows = sum(ageStateChange(end - 1, 2:end-1));
    
    if(birthEventsFromNonInfectedCows + birthEventsFromInfectedCows > 40)
        fprintf("[WARN]: More than 40 birth events are being attempted in this call of" + ...
        " Farm_Square_Model.Update_Cattle_State(). This may indicate that τ is too large.\n");
    end
    
    if femaleBirthsOnly

        % Draw the numbers of successful births from non-infected/infected cows from the
        % binomial distribution with success probability par.eta * 0.5.
        successfulBirthsFromNonInfectedCows = binornd(birthEventsFromNonInfectedCows, ...
        par.eta * 0.5);
        successfulBirthsFromInfectedCows = binornd(birthEventsFromInfectedCows, ...
        par.eta * 0.5);
    else

        % Draw the numbers of successful births from non-infected/infected cows from the
        % binomial distribution with success probability par.eta.
        successfulBirthsFromNonInfectedCows = binornd(birthEventsFromNonInfectedCows, ...
        par.eta);
        successfulBirthsFromInfectedCows = binornd(birthEventsFromInfectedCows, ...
        par.eta);
    end
    
    % Draw the number of the successful births from infected cows in which transplacental
    % transmission occurs from the binomial distribution with success probability par.zeta.
    tptBirthsFromInfectiousCows = binornd(successfulBirthsFromInfectedCows, par.zeta);
    
    % Add non-infected successful births into the X_Y1 compartment.
    cattleStateMatrix(1,1) = cattleStateMatrix(1,1) + ...
    successfulBirthsFromNonInfectedCows + (successfulBirthsFromInfectedCows - ...
    tptBirthsFromInfectiousCows);
    
    % Add infected successful births into the Y_1,Y1 compartment.
    cattleStateMatrix(1,2) = cattleStateMatrix(1,2) + tptBirthsFromInfectiousCows;
    
    % ----- END SECTION -----
end