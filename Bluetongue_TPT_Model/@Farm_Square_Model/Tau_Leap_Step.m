function exitCode = Tau_Leap_Step(self, ratesRequireRecalculation, tau, temp)

    % Evolves the Farm_Square_Model object according in time according to the "tau-leap"
    % algorithm. The function takes a single tau-leap step and updates
    % self.dairyCattleStates, self.beefCattleStates, self.sheepStates and
    % self.vectorStates based on Poisson-distributed draws to determine the number of
    % different types of events (infection, aging, conception, birth, culling, etc.) which
    % are assumed to occur in the time step.
    %
    % ratesRequireRecalculation: A vector containing 6 boolean values. Each value
    % specifies whether a set of event rate matrices should be calculated this time step.
    % If false, a previously-calculated "cached" version of the rate matrix is used instead.
    % The 6 values correspond to the following matrices (in order):
    % Bool 1: dairyCattleInfectionRateMatrix, beefCattleInfectionRateMatrix,
    % sheepInfectionRateMatrix.
    % Bool 2: dairyCattleAgeRateMatrix, beefCattleAgeRateMatrix.
    % Bool 3: dairyCattleMortalityRateMatrix, beefCattleMortalityRateMatrix,
    % sheepMortalityRateMatrix.
    % Bool 4: dairyCattleCullingRateMatrix, beefCattleCullingRateMatrix,
    % sheepCullingRateMatrix.
    % Bool 5: vectorInfectionRateMatrix.
    % Bool 6: vectorMortalityRateMatrix.
    %
    % tau: The time interval (in days) to tau-leap by.
    % temp: The current farm square temperature.
    %
    % Returns an 'exit code', which is one of the following integers:
    % 0: Tau-leap step taken successfully and infection is present in the farm square.
    % 1: Tau-leap step taken successfully but infection is not present in the farm square.
    % 2: There are no cattle or sheep in the farm square. Tau-leap step not taken.
    %
    % AUTHOR: Laurence Dhonau.

    par = self.parameters;

    % Initialise an integer to store the number of times we have drawn a set of event
    % counts from the Poisson distributions.
    drawCount = 0;

    % Initialise a bool to indicate whether we have successfully drawn event counts yet
    % (i.e. drawn counts such that the drawn counts would not make any of the population
    % classes become negative) and subsequently updated the model state.
    modelStateUpdated = false;
 
    % Draw event counts until we successfully update the model state.
    while (~modelStateUpdated)
     
        % Increment the number of Poisson draws by 1.
        drawCount = drawCount + 1;

        % Compute the current total dairy cattle, beef cattle, sheep and vector
        % populations.
        H_D = sum(sum(self.dairyCattleStates.tensor(:,:,end)));
        H_B = sum(sum(self.beefCattleStates.tensor(:,:,end)));
        H_S = sum(self.sheepStates.tensor(:,:,end));
        N = sum(self.vectorStates.tensor(:,:,end));

        % If there are no cattle or sheep in the farm square then return immediately with
        % exit code 2.
        if sum([H_D H_B H_S]) == 0
            fprintf("[WARN] Tau_Leap_Step was called with no hosts in the farm square.\n");
            exitCode = 2;
            return;
        end

        % Compute the current proportion of vector bites on cattle.
        phi = (H_D + H_B) / (H_D + H_B + par.sigma * H_S);

        % Compute the current infected cattle, sheep and vector totals.
        infectedCattleCount = sum(sum(self.dairyCattleStates.tensor(:,2:end-1,end))) + ...
        sum(sum(self.beefCattleStates.tensor(:,2:end-1,end)));
        infectedSheepCount = sum(self.sheepStates.tensor(1, 2:end-1, end));
        infectedVectorCount = self.vectorStates.tensor(1, end, end);

        %% ----- SECTION: EVENT RATES -----

        if ratesRequireRecalculation(1)

            % Compute the host infection transition matrices.
            [self.dairyCattleInfectionRateMatrix, self.beefCattleInfectionRateMatrix, ...
            self.sheepInfectionRateMatrix] = self.Compute_Host_Infection_Rates(par, ...
            self.t, temp, H_D, H_B, H_S, infectedVectorCount / N);
        end

        if ratesRequireRecalculation(2)

            % Compute the cattle age transition matrices.
            [self.dairyCattleAgeRateMatrix, self.beefCattleAgeRateMatrix] = ...
            self.Compute_Cattle_Age_Rates(par);
        end

        if ratesRequireRecalculation(3)

            % Compute the host mortality matrices.
            [self.dairyCattleMortalityRateMatrix, self.beefCattleMortalityRateMatrix, ...
            self.sheepMortalityRateMatrix] = self.Compute_Host_Mortality_Rates(par);
        end

        if ratesRequireRecalculation(4)

            % Compute the host culling matrices.
            [self.dairyCattleCullingRateMatrix, self.beefCattleCullingRateMatrix, ...
            self.sheepCullingRateMatrix] = self.Compute_Host_Culling_Rates(par);
        end

        if ratesRequireRecalculation(5)

            % Compute the vector infection transition matrix.
            if H_D + H_B == 0
                self.vectorInfectionRateMatrix = self.Compute_Vector_Infection_Rates(par, ...
                self.t, temp, phi, 0, infectedSheepCount / H_S);
            elseif H_S == 0
                self.vectorInfectionRateMatrix = self.Compute_Vector_Infection_Rates(par, ...
                self.t, temp, phi, infectedCattleCount / (H_D + H_B), 0);
            else
                self.vectorInfectionRateMatrix = self.Compute_Vector_Infection_Rates(par, ...
                self.t, temp, phi, infectedCattleCount / (H_D + H_B), infectedSheepCount / H_S);
            end
        end

        if ratesRequireRecalculation(6)

            % Compute the vector mortality matrix.
            self.vectorMortalityRateMatrix = self.Compute_Vector_Mortality_Rates(par, temp);
        end

        % ----- END SECTION -----

        %% ----- SECTION: EVENT DRAWS -----

        % Draw the event counts for dairy cattle events in this time step from the Poisson
        % distribution.
        dairyCattleInfectionStateChange = poissrnd(tau * self.dairyCattleInfectionRateMatrix .* self.dairyCattleStates.tensor(:,:,end));
        dairyCattleAgeStateChange = poissrnd(tau * self.dairyCattleAgeRateMatrix .* self.dairyCattleStates.tensor(:,:,end));
        dairyCattleMortalityStateChange = poissrnd(tau * self.dairyCattleMortalityRateMatrix .* self.dairyCattleStates.tensor(:,:,end));
        dairyCattleCullingStateChange = poissrnd(tau * self.dairyCattleCullingRateMatrix .* self.dairyCattleStates.tensor(:,:,end));
        
        % Draw the event counts for beef cattle events in this time step from the Poisson
        % distribution.
        beefCattleInfectionStateChange = poissrnd(tau * self.beefCattleInfectionRateMatrix .* self.beefCattleStates.tensor(:,:,end));
        beefCattleAgeStateChange = poissrnd(tau * self.beefCattleAgeRateMatrix .* self.beefCattleStates.tensor(:,:,end));
        beefCattleMortalityStateChange = poissrnd(tau * self.beefCattleMortalityRateMatrix .* self.beefCattleStates.tensor(:,:,end));
        beefCattleCullingStateChange = poissrnd(tau * self.beefCattleCullingRateMatrix .* self.beefCattleStates.tensor(:,:,end));

        % Draw the event counts for sheep events in this time step from the Poisson
        % distribution.
        sheepBirthStateChange = poissrnd(tau * par.b_S * H_S);
        sheepInfectionStateChange = poissrnd(tau * self.sheepInfectionRateMatrix .* self.sheepStates.tensor(:,:,end));
        sheepMortalityStateChange = poissrnd(tau * self.sheepMortalityRateMatrix .* self.sheepStates.tensor(:,:,end));
        sheepCullingStateChange = poissrnd(tau * self.sheepCullingRateMatrix .* self.sheepStates.tensor(:,:,end));

        % Draw the event counts for vector events in this time step from the Poisson
        % distribution.
        vectorInfectionStateChange = poissrnd(tau * self.vectorInfectionRateMatrix .* self.vectorStates.tensor(:,:,end));
        vectorMortalityStateChange = poissrnd(tau * self.vectorMortalityRateMatrix .* self.vectorStates.tensor(:,:,end));

        % ----- END SECTION -----

        %% ----- SECTION: UPDATE CATTLE STATES -----

        % Provisionally update the dairy cattle state matrix based on the drawn event
        % counts, handling any birth events (including potential transplacental
        % transmission), appending the new matrix to the end of the dairyCattleStates
        % tensor.
        self.dairyCattleStates.tensor(:,:,end + 1) = self.Update_Cattle_State(par.L, ...
        par.M, self.dairyCattleStates.tensor(:,:,end), dairyCattleInfectionStateChange, ...
        dairyCattleAgeStateChange, dairyCattleMortalityStateChange, ...
        dairyCattleCullingStateChange, self.dairyAgeReEntryRow, true);

        % Provisionally update the beef cattle state matrix based on the drawn event
        % counts, handling any birth events (including potential transplacental
        % transmission), appending the new matrix to the end of the beefCattleStates
        % tensor.
        self.beefCattleStates.tensor(:,:,end + 1) = self.Update_Cattle_State(par.Ltilde, ...
        par.M, self.beefCattleStates.tensor(:,:,end), beefCattleInfectionStateChange, ...
        beefCattleAgeStateChange, beefCattleMortalityStateChange, ...
        beefCattleCullingStateChange, self.beefAgeReEntryRow, false);

        % ----- END SECTION -----

        %% ----- SECTION: UPDATE VECTOR STATES -----

        % Provisionally update the vector state matrix based on drawn infection event
        % counts.
        self.vectorStates.tensor(:,:,end + 1) = self.vectorStates.tensor(:, :, end) + ...
        [0 vectorInfectionStateChange(1:end-1)] - vectorInfectionStateChange;

        % Further update the vector state matrix based on drawn mortality event counts.
        self.vectorStates.tensor(:,:,end) = self.vectorStates.tensor(:, :, end) - ...
        vectorMortalityStateChange;

        % We model the total vector population as constant so we must include compensatory
        % recruitment to balance any loss in vectors (mortality). The newly-recruited
        % vectors enter the suspectible compartment.
        self.vectorStates.tensor(:,1,end) = self.vectorStates.tensor(:, 1, end) + ...
        sum(vectorMortalityStateChange);

        % ----- END SECTION -----

        %% ----- SECTION: UPDATE SHEEP STATES -----

        % Provisionally update the sheep state matrix based on drawn infection event
        % counts.
        self.sheepStates.tensor(:,:,end + 1) = self.sheepStates.tensor(:, :, end) + ...
        [0 sheepInfectionStateChange(1:end-1)] - sheepInfectionStateChange;

        % Further update the sheep state matrix based on drawn mortality and culling event
        % counts.
        self.sheepStates.tensor(:,:,end) = self.sheepStates.tensor(:, :, end) - ...
        sheepMortalityStateChange - sheepCullingStateChange;

        % Further update the sheep state matrix based on the drawn birth event count.
        self.sheepStates.tensor(:,1,end) = self.sheepStates.tensor(:,1,end) + ...
        sheepBirthStateChange;

        % ----- END SECTION -----

        %% ----- SECTION: VERIFY MODEL STATE UPDATE -----

        if max(max(isnan(self.dairyCattleStates.tensor(:,:,end)))) == 1 ... 
        || max(max(isnan(self.beefCattleStates.tensor(:,:,end)))) == 1 ...
        || max(isnan(self.sheepStates.tensor(:,:,end))) == 1 ...
        || max(isnan(self.vectorStates.tensor(:,:,end))) == 1
            input("[ERROR] Tau_Leap_Step has produced at least one NaN value for the new" + ...
            " model state. Something has gone very wrong. Press Return to continue" + ...
            " (or use a breakpoint to pause here). ");
        end

        % Check if the provisionally-updated model state has all non-negative compartment
        % values.
        if (min(min(self.dairyCattleStates.tensor(:,:,end))) >= 0) ... 
        && (min(min(self.beefCattleStates.tensor(:,:,end))) >= 0) ...
        && (min(self.vectorStates.tensor(:,:,end)) >= 0) ...
        && (min(self.sheepStates.tensor(:,:,end)) >= 0)
         
            % We have successfully updated the model state so we can exit the loop.
            modelStateUpdated = true;
        else

            % Check if we have already drawn 2 sets of event counts.
            %if drawCount == 2

                % Warn the user that we will now perform a 3rd set of Poisson draws. If
                % this happens for many time steps, tau-leap will be very slow and tau is
                % likely too large.
                %fprintf("[WARN] Trying τ-leap Poisson draws for the 3rd time. This" + ...
                %" may indicate τ is too large.\n");
            %end

            % Check if we have already drawn 15 sets of event counts.
            if drawCount == 15
                input(sprintf("[ERROR] τ-leap is likely stuck at t = %0.1f days (have tried 15" + ...
                " Poisson draws). This may indicate τ is too large. Press Return to continue" + ...
                " (or use a breakpoint to pause here). ", self.t));
            end

            % Discard the model state updates we attempted.
            self.dairyCattleStates.tensor(:, :, end) = [];
            self.beefCattleStates.tensor(:, :, end) = [];
            self.sheepStates.tensor(:, :, end) = [];
            self.vectorStates.tensor(:, :, end) = [];
        end
    end

    % ----- END SECTION -----
 
    % Update the Farm_Square_Model object time
    self.t = self.t + tau;

    % Compute the new infected cattle, sheep and vector totals.
    infectedCattleCount = sum(sum(self.dairyCattleStates.tensor(:,2:end-1,end))) + ...
    sum(sum(self.beefCattleStates.tensor(:,2:end-1,end)));
    infectedSheepCount = sum(self.sheepStates.tensor(1, 2:end-1, end));
    infectedVectorCount = self.vectorStates.tensor(1, end, end);

    % If there is no infection in the farm square.
    if sum([infectedCattleCount, infectedSheepCount, infectedVectorCount]) == 0
        exitCode = 1;
    else
        exitCode = 0;
    end
end