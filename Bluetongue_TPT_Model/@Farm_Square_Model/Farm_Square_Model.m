classdef Farm_Square_Model < handle

    % Encapsulates information about a single farm square and its host/vector populations.
    % Also contains the properties and methods required to simulate BTV transmission
    % within the square.
    %
    % AUTHOR: Laurence Dhonau.

    properties
        parameters;
        farmLatitude;
        farmLongitude;
        initialPopulations;

        dairyCattleStates;
        beefCattleStates;
        sheepStates;
        vectorStates;

        t;
        dairyAgeReEntryRow;
        beefAgeReEntryRow;

        dairyCattleInfectionRateMatrix;
        beefCattleInfectionRateMatrix;
        sheepInfectionRateMatrix;
        dairyCattleAgeRateMatrix;
        beefCattleAgeRateMatrix;
        dairyCattleMortalityRateMatrix;
        beefCattleMortalityRateMatrix;
        sheepMortalityRateMatrix;
        dairyCattleCullingRateMatrix;
        beefCattleCullingRateMatrix;
        sheepCullingRateMatrix;
        vectorInfectionRateMatrix;
        vectorMortalityRateMatrix;
    end

    methods(Static)

        %% ----- SECTION: STATIC METHOD SIGNATURES -----

        % @Farm_Square_Model/Build_Cattle_Label_Matrix.m
        cattleLabelMatrix = Build_Cattle_Label_Matrix(classCounts, classLabels);

        % @Farm_Square_Model/Compute_Cattle_Age_Rates.m
        [dairyCattleAgeRateMatrix, beefCattleAgeRateMatrix] = Compute_Cattle_Age_Rates(par);

        % @Farm_Square_Model/Compute_Host_Culling_Rates.m
        [dairyCattleCullingRateMatrix, beefCattleCullingRateMatrix, ...
        sheepCullingRateMatrix] = Compute_Host_Culling_Rates(par);

        % @Farm_Square_Model/Compute_Host_Infection_Rates.m
        [dairyCattleInfectionRateMatrix, beefCattleInfectionRateMatrix, ...
        sheepInfectionRateMatrix] = Compute_Host_Infection_Rates(par, t, temp, H_D, H_B, ...
        H_S, infectedVectorFraction);

        % @Farm_Square_Model/Compute_Host_Mortality_Rates.m
        [dairyCattleMortalityRateMatrix, beefCattleMortalityRateMatrix, ...
        sheepMortalityRateMatrix] = Compute_Host_Mortality_Rates(par);

        % @Farm_Square_Model/Compute_Vector_Activity.m
        activity = Compute_Vector_Activity(par, t, activityNormaliser)

        % @Farm_Square_Model/Compute_Vector_Infection_Rates.m
        vectorInfectionRateMatrix = Compute_Vector_Infection_Rates(par, t, temp, phi, ...
        infectedCattleFraction, infectedSheepFraction)

        % @Farm_Square_Model/Compute_Vector_Mortality_Rates.m
        vectorMortalityRateMatrix = Compute_Vector_Mortality_Rates(par, temp)
    
        % ----- END SECTION -----
    end

    methods
        %% ----- SECTION: OBJECT CONSTRUCTOR -----

        function self = Farm_Square_Model(farmLatitude, farmLongitude, parameters, ...
        initialTime, initialPopulations)

            % Constructs a Farm_Square_Model object with the passed values. In particular,
            % the function sets the dairyCattleStates, beefCattleStates, sheepStates and
            % vectorStates properties to Labelled_Tensor objects with the initial
            % conditions specified in initialPopulations (initial conditions can be
            % specified in various ways; see Farm_Square_Model.Define_Initial_Conditions()).
            %
            % farmLatitude: The latitude of the centre of the farm square (in decimal
            % degrees) to 3 decimal places.
            % farmLongitude: The longitude of the centre of the farm square (in decimal
            % degrees) to 3 decimal places.
            % parameters: A structure specifying within-farm transmission parameters for
            % use in the model. Bluetongue_TPT_Model.Define_Within_Farm_Parameters()
            % generates a structure of the required form.
            % initialTime: The value of time (t), in days, to initialise the model with.
            % initialPopulations: A struct with three fields ("dairy", "beef", "sheep")
            % specifying the host population structures to initialise the model with. The
            % full form of this structure is specified in Farm_Square_Model.Define_
            % Initial_Conditions().
            
            generatedValidParams = false;

            % Draw vector-to-host ratio parameters until they are positive. Note that the
            % probability we have to redraw is low (< 10%).
            while ~generatedValidParams

                % Compute vector-to-host ratio for cattle in the farm square.
                parameters.m_C = gamrnd(parameters.s_V, parameters.mu_V / parameters.s_V);
            
                % Compute vector-to-host ratio for sheep in the farm square.
                parameters.m_S = gamrnd(parameters.s_V, parameters.mu_V / parameters.s_V);

                if parameters.m_C >= 0 && parameters.m_S >= 0
                    generatedValidParams = true;
                end
            end

            self.farmLatitude = farmLatitude;
            self.farmLongitude = farmLongitude;
            self.parameters = parameters;
            self.t = initialTime;
            self.initialPopulations = initialPopulations;

            % Instantiate a Labelled_Tensor to store the evolution of the dairy cattle
            % state matrix, setting up the label matrix to accurately reflect the
            % compartments described in the project report and setting the tensor to a
            % zero matrix initially. Store the result in the dairyCattleStates property.
            dairyCattleStatesLabelMatrix = Farm_Square_Model.Build_Cattle_Label_Matrix ...
            ([self.parameters.L, self.parameters.M, self.parameters.K], ["Y", "A", "P", "V"]);
            dairyCattleStatesTensor = zeros([size(dairyCattleStatesLabelMatrix), 1]);
            self.dairyCattleStates = Labelled_Tensor(dairyCattleStatesTensor, ...
            dairyCattleStatesLabelMatrix);

            % Instantiate a Labelled_Tensor for beef cattle similarly.
            beefCattleStatesLabelMatrix = Farm_Square_Model.Build_Cattle_Label_Matrix ...
            ([self.parameters.Ltilde, self.parameters.M, self.parameters.K], ["Y", "A", "P", "V"]);
            beefCattleStatesTensor = zeros([size(beefCattleStatesLabelMatrix), 1]);
            self.beefCattleStates = Labelled_Tensor(beefCattleStatesTensor, ...
            beefCattleStatesLabelMatrix);

            % Instantiate a Labelled_Tensor for sheep similarly.
            sheepStatesLabelMatrix = ["X.", vecsprintf("Y_%s.", string(split(num2str(1 : ...
            self.parameters.Kbar)))).', "Z."];
            sheepStatesTensor = zeros([size(sheepStatesLabelMatrix), 1]);
            self.sheepStates = Labelled_Tensor(sheepStatesTensor, sheepStatesLabelMatrix);

            % Instantiate a Labelled_Tensor for vectors similarly.
            vectorStatesLabelMatrix = ["S.", vecsprintf("L_%s.", string(split(num2str(1 : ...
            self.parameters.k)))).', "I."];
            vectorStatesTensor = zeros([size(vectorStatesLabelMatrix), 1]);
            self.vectorStates = Labelled_Tensor(vectorStatesTensor, vectorStatesLabelMatrix);

            % Apply the initial conditions specified in self.initialPopulations to the
            % host/vector state Labelled_Tensors we have just instantiated.
            self.Define_Initial_Conditions();

            % Set the re-entry row properties so dairy/beef cattle exit postpartum classes
            % into the corresponding adult class (the (L + 1) th row in the state
            % matrices), if they are not culled.
            self.dairyAgeReEntryRow = self.parameters.L + 1;
            self.beefAgeReEntryRow = self.parameters.Ltilde + 1;
        end

        % ----- END SECTION -----

        %% ----- SECTION: METHOD SIGNATURES -----

        % @Farm_Square_Model/Define_Initial_Conditions.m
        Define_Initial_Conditions(self);

        % @Farm_Square_Model/Tau_Leap_Step.m
        exitCode = Tau_Leap_Step(self, ratesRequireRecalculation, tau, temp);

        % @Farm_Square_Model/Update_Cattle_State.m
        cattleStateMatrix = Update_Cattle_State(self, L, M, cattleStateMatrix, ...
        infectionStateChange, ageStateChange, mortalityStateChange, cullingStateChange, ...
        ageReEntryRow, femaleBirthsOnly);

        % ----- END SECTION -----
    end
end