function Run(TPTVALUE)

    % Run the BTV model across all of France.
    %
    % AUTHOR: Laurence Dhonau.

    % We initialise the farmSquareDistanceMatrix as a global variable since it is very
    % large (approx. 600 MB) so takes a non-negligible amount of time to load into
    % memory. By using a global variable, we only need load the matrix into memory
    % on the first simulation after launching MATLAB.
    global farmSquareDistanceMatrix;

    % Generate a name used for the various export files.
    exportFileName = sprintf("BTV_Simulation_%s_%d", datetime('today'), ...
    round(rand * 100000));

    % Create an Outputs folder if it does not already exist.
    [~,~] = mkdir("Outputs");

    % Save all Command Window output to a log file until the function exits.
    diary(sprintf("Outputs/%s.log", exportFileName));

    fprintf("[INFO] This is the Bluetongue_TPT_Model v1.0. Starting simulation.\n");

    %% ----- SECTION: SIMULATION PARAMETERS -----

    % Define the start and end times for the simulation (inclusive). Time is measured in
    % days, where day 0 is 01/01/2006. The end time must not exceed day 1460 (31/12/2009)
    % since we do not import temperature data for days after this time.
    simulationStartTime = 236;
    simulationEndTime = 236 + 300;

    % Define the time step for the simulation. This must currently be set to 1 (otherwise
    % expect many errors) due to the implementation of between-farm diffusion we inherit
    % from the base model.
    tau = 1;

    % Define the scale factor for within-farm tau-leap. The value of "tau" used for
    % tau-leap is given by tau * tauScale. Note that tauScale should be the reciprocal of
    % a positive integer.
    tauScale = 1/20;

    % Define the farm squares used to seed BTV transmission. Each row in
    % seedFarmSquareLongLats should contain two elements: the longitude and latitude (in
    % decimal degrees to 3 decimal places) of a farm square. Each farm square used to seed
    % the simulation must contain at least 1 host in the GLW v3 dataset.
    seedFarmSquareLongLats = [4.708, 49.958;
                              4.792, 50.042];

    % Define the parameters used for between-farm BTV transmission: gamma (transmission
    % parameter) and D (diffusion coefficient).
    parameters = struct("gamma", betarnd(3.49, 2.14), "D", normrnd(3.01, 1.01));

    % Define the bounds for the Uniform distribution from which the fraction of cattle we
    % assume to be dairy cattle (rather than beef cattle) in each farm square is drawn.
    dairyCattleFractionBounds = [0.38, 0.52];

    % Define the bounds for the Uniform distribution from which the number of infected
    % vectors we insert on a farm square when it becomes infected is drawn.
    initialInfectedVectorCountBounds = [1, 100];

    % Define the size of the MATLAB parallel pool to use for tau-leap.
    parallelPoolSize = 16;

    fprintf("[INFO] Drawing within-farm transmission parameters.\n");

    generatedValidParams = false;

    % Draw within-farm transmission parameters until the vector activity parameters are
    % non-positive and all other parameters are non-negative. Note that the probability we
    % have to redraw is low (< 10%).
    while ~generatedValidParams

        % Construct a structure containing within-farm transmission parameters.
        farmSquareParameters = Bluetongue_TPT_Model.Define_Within_Farm_Parameters();

        % If all parameters are of valid sign.
        if min(cell2mat(struct2cell(rmfield(farmSquareParameters, {'b_11', 'b_12', 'b_21', ...
        'b_22'})))) >= 0 && max([farmSquareParameters.b_11, farmSquareParameters.b_12, ...
        farmSquareParameters.b_21, farmSquareParameters.b_22]) <= 0
            generatedValidParams = true;
        else
            fprintf("[WARN] Invalid within-farm parameters drawn. Redrawing...\n");
        end
    end

    % ----- END SECTION -----

    %% ----- SECTION: SIMULATION DATASETS -----

    fprintf("[INFO] Loading required datasets.\n");

    % Load a table containing farm square information. See Bluetongue_TPT_Model.Construct_
    % Farm_Square_Table() for construction details.
    farmSquareTable = Bluetongue_TPT_Model.Load_Dataset("Datasets/Farm_Square_Table.mat");

    % Only load farmSquareDistanceMatrix from disk if it is not already in memory (from a
    % previous simulation).
    if max(size(farmSquareDistanceMatrix)) == 0

        % Load a matrix containing the pairwise distance between each of the 9122 farm
        % squares. See Bluetongue_TPT_Model.Compute_Lat_Long_Distance_Matrix() for
        % construction details.
        farmSquareDistanceMatrix = Bluetongue_TPT_Model.Load_Dataset(...
        "Datasets/Farm_Square_Pairwise_Dist.mat");
    end

    % Load a dictionary containing information on which weather station is closest in
    % distance to each farm square.
    farmSquareStationDictionary = Bluetongue_TPT_Model.Load_Dataset(...
    "Datasets/Farm_Square_Nearest_Station_Dictionary.mat");

    % Load a dictionary containing the minimum/maximum daily temperatures for each
    % weather station.
    stationTemperatureDictionary = Bluetongue_TPT_Model.Load_Dataset(...
    "Datasets/Station_Temperature_Dictionary.mat");

    % Load a dictionary containing giving the vector mortality rate for each 0.1°C
    % temperature increment in the range -50 °C to 50 °C inclusive. See Bluetongue_TPT_
    % Model.Construct_Vector_Mortality_Dictionary() for construction details.
    vectorMortalityDictionary = Bluetongue_TPT_Model.Load_Dataset(...
    "Datasets/Vector_Mortality_Dictionary.mat");

    % ----- END SECTION -----

    %% ----- SECTION: SIMULATION SEEDING -----

    fprintf("[INFO] Seeding BTV simulation.\n");

    % Compute the number of farm squares.
    farmSquareCount = size(farmSquareTable, 1);

    % Initalise a cell array to keep track of farm squares which have been infected at
    % some point during the simulation (but may not be infected currently). The cell array
    % consists of cells containing two elements: the ID of the farm square and a
    % Farm_Square_Model object.
    infectedFarmSquareModels = {};

    % Iterate over the farm squares we will seed the simulation with.
    for i = 1:size(seedFarmSquareLongLats, 1)

        % Retrieve the ID of the farm square to be seeded based on its latitude and
        % longitude.
        farmSquareID = Bluetongue_TPT_Model.Lookup_Farm_Square_Lat_Long(farmSquareTable, ...
        seedFarmSquareLongLats(i,2),seedFarmSquareLongLats(i,1));

        % Retrieve the row for the farm square from the farm square table.
        farmSquareRow = table2array(farmSquareTable(farmSquareID, :));

        % Draw the fraction of cattle in the farm square which are dairy cattle from the
        % Uniform distribution.
        dairyCattleFraction = rand * (dairyCattleFractionBounds(2) - ...
        dairyCattleFractionBounds(1)) + dairyCattleFractionBounds(1);

        % Instantiate a new Farm_Square_Model object and set it to a demographic
        % steady state with host populations from the farm square table.
        farmSquareModel = Bluetongue_TPT_Model.Initialise_Farm_Square_Model( ...
        farmSquareRow(2), farmSquareRow(3), farmSquareParameters, simulationStartTime, ...
        round([dairyCattleFraction * farmSquareRow(4), (1 - dairyCattleFraction) * ...
        farmSquareRow(4), farmSquareRow(5)]));

        % Draw the number of infected vectors to insert into the farm square from the
        % Uniform distribution.
        initialInfectedVectorCount = round(rand * (initialInfectedVectorCountBounds(2) - ...
        initialInfectedVectorCountBounds(1)) + initialInfectedVectorCountBounds(1));

        % Insert infected vectors to seed infection *within* the farm square.
        farmSquareModel.vectorStates.Set_Val("I.", 1, initialInfectedVectorCount);

        % Update the farm square table to indicate that the farm square has been
        % infected, along with the time of infection.
        farmSquareTable(farmSquareID, 6) = {true};
        farmSquareTable(farmSquareID, 7) = {simulationStartTime};
        farmSquareTable(farmSquareID, 8) = {simulationStartTime};

        % Add the farm square ID and model object to the infected squares.
        infectedFarmSquareModels{end + 1} = {farmSquareID, farmSquareModel};
    end

    % ----- END SECTION -----

    %% ----- SECTION: BETWEEN-FARM TRANSMISSION AND TAU-LEAPING -----

    fprintf("[INFO] Simulation seeded! Initialising a parallel pool.\n");

    % Initialise a parallel pool with the requested number of workers
    if isempty(gcp('nocreate'))
        parpool('Threads', parallelPoolSize);
    end

    fprintf("[INFO] Parallel pool initialised. Beginning simulation.\n");

    simulationTime = simulationStartTime;

    % Initalise a vector to store a time series of the number of infected farms at each
    % point in time.
    infectedFarmTimeSeries = [size(infectedFarmSquareModels, 2) zeros(1, simulationEndTime - simulationStartTime)];

    % Construct a dictionary where the keys are each time (t) value featured in
    % the simulation and the values are the vector activity values on each day (t
    % value).
    vectorActivityDictionary = dictionary(simulationStartTime : simulationEndTime, ...
    Farm_Square_Model.Compute_Vector_Activity(farmSquareParameters, simulationStartTime : ...
    simulationEndTime, farmSquareParameters.activityNormaliser));

    % Convert the farm square table into a matrix while performing tau-leap since this
    % *dramatically* reduces overall computational cost.
    farmSquareTableMatrix = table2array(farmSquareTable);

    % Start a stopwatch to keep track of the time spent running the simulation.
    simulationStopwatch = tic;

    while(simulationTime <= simulationEndTime)

        fprintf("[INFO] ----- Running simulation for day %d. -----\n", simulationTime);

        % Allow the user to pause the simulation if they wish every 10 time steps.
        if rem(simulationTime - simulationStartTime, 10) == 0 && isfile("PAUSE")
            input("[INPUT] Simulation paused. Press Return to continue. ");
        end

        %% ----- SUBSECTION: BETWEEN-FARM TRANSMISSION -----

        betweenFarmTransmissionStopwatch = tic;

        % Initialise a 'cache' cell array to store key information about elements of
        % infectedFarmSquareModels which we will access repeatedly in this time step while
        % determining whether each infected farm infects each non-infected one.
        infectedFarmSquareCache = cell(size(infectedFarmSquareModels));
        
        % Iterate over the elements of infectedFarmSquareModels.
        for j = 1:size(infectedFarmSquareModels, 2)
            infectedFarmSquareID = infectedFarmSquareModels{j}{1};

            % Retrieve a bool indicating whether the farm square is currently infected.
            farmSquareInfected = farmSquareTableMatrix(infectedFarmSquareID, 6);

            if farmSquareInfected

                % Retrieve the time of most recent infection for the farm square.
                tau_j = farmSquareTableMatrix(infectedFarmSquareID, 8);

                % Retrieve the time of first infection for the farm square.
                taudash_j = farmSquareTableMatrix(infectedFarmSquareID, 7);

                % Compute a time series for the midpoint of the minimum and maximum daily
                % vector mortalities for the farm square between tau_j and simulationTime
                % (inclusive). We will use this repeatedly this time step.
                farmSquareMidMortalityTimeSeries = Bluetongue_TPT_Model. ...
                Compute_Mid_Mortality_Time_Series(infectedFarmSquareID, tau_j, ...
                simulationTime, farmSquareStationDictionary, stationTemperatureDictionary, ...
                vectorMortalityDictionary);

                % Retrieve the row of the farmSquareDistanceMatrix corresponding to the
                % farm square.
                farmSquareDistances = farmSquareDistanceMatrix(infectedFarmSquareID, :);
    
                % Add the mortality time series to the cache, along with the square ID,
                % row of farmSquareDistanceMatrix and times of most recent and first
                % infection.
                infectedFarmSquareCache{j} = {farmSquareInfected, infectedFarmSquareID, ...
                tau_j, farmSquareMidMortalityTimeSeries, taudash_j, farmSquareDistances};
            else

                % If the farm square is not currently infected, we need not cache
                % information and we add false to the cache to indicate this.
                infectedFarmSquareCache{j} = {farmSquareInfected};
            end
        end

        % Store the size of infectedFarmSquareModels since this will increase if farm
        % squares are infected for the first time during the time step.
        sizeOfInfectedFarmSquareModels = size(infectedFarmSquareModels, 2);

        % Initialise a logical vector to track whether each farm square gains infection in
        % this time step.
        squareGainedInfection = false(1, farmSquareCount);

        % Iterate over all farm squares using parallel workers.
        parfor k = 1:farmSquareCount

            % Retrieve the row for the farm square from the farm square table.
            farmSquareRow = farmSquareTableMatrix(k, :);

            % If the farm square has at least 1 host and is not currently infected.
            if ~farmSquareRow(6) && sum(farmSquareRow(4:5)) > 0 && ~farmSquareRow(9)

                % Iterate over all the potentially infected farm squares.
                for j = 1:sizeOfInfectedFarmSquareModels

                    % If the square is not currently infected then skip it.
                    if ~infectedFarmSquareCache{j}{1}
                        continue;
                    end

                    % Get the ID and Farm_Square_Model object for the potentially infected
                    % square.
                    infectedFarmSquareID = infectedFarmSquareModels{j}{1};
                    infectedFarmSquareModel = infectedFarmSquareModels{j}{2};

                    % Before making use of the infected farm square cache, check that the
                    % square ID of the jth element matches up with the ID from
                    % infectedFarmSquareModels (for safety).
                    if infectedFarmSquareCache{j}{2} ~= infectedFarmSquareID
                        input("[ERROR] Mismatch between infectedFarmSquareCache and" + ...
                        " infectedFarmSquareModels. Something has gone very wrong. Press" + ...
                        " Return to continue (or use a breakpoint to pause here). ");
                    end

                    % Retrieve the distance between the non-infected farm square and the
                    % infected one.
                    farmSquareSeparation = infectedFarmSquareCache{j}{6}(k);

                    % Retrieve the t value on which the infected square became infected
                    % (most recently).
                    tau_j = infectedFarmSquareCache{j}{3};

                    % Retrieve the t value on which the infected square was first infected.
                    taudash_j = infectedFarmSquareCache{j}{5};

                    if exp(-farmSquareSeparation ^ 2 / (4 * parameters.D * (simulationTime ...
                    - tau_j + 1))) == 0

                        % In this case, the farm square separation is too great, or the
                        % farm square infection too recent, for infection to have spread
                        % from the infected farm square to the non-infected one so we
                        % continue to the next infected square.
                        continue;
                    end

                    farmSquareMidMortalityTimeSeries = infectedFarmSquareCache{j}{4};

                    rateSum = 0;

                    % Iterate over the days since the farm square became infected.
                    for tdash = tau_j:simulationTime

                        % Sum up the BTV diffusion terms (using the form given in the base
                        % model).
                        rateSum = rateSum + vectorActivityDictionary(tdash) * ...
                        infectedFarmSquareModel.vectorStates.tensor(1, end, tdash - taudash_j ...
                        + 1) / (4 * pi * parameters.D * (simulationTime - tdash + 1)) * ...
                        exp(-farmSquareSeparation ^ 2 / (4 * parameters.D * ...
                        (simulationTime - tdash + 1))) * exp(-sum( ...
                        farmSquareMidMortalityTimeSeries(tdash - tau_j + 1 : simulationTime - ...
                        tau_j + 1)));
                    end

                    % Compute the probability that the infected farm infects the
                    % non-infected one in this time step using the sum of the diffusion
                    % terms.
                    probabilityOfInfection = 1 - exp(-parameters.gamma * rateSum);

                    % Perform a Bernoulli trial to determine if transmission occurs.
                    if rand <= probabilityOfInfection

                        % Transmission does occur so update squareGainedInfection and
                        % proceed to the next non-infected farm square.
                        squareGainedInfection(k) = true;
                        break;
                    end
                end
            end
        end

        % Retrieve a vector of the farm squares which gained infection this time step.
        newlyInfectedSquareIDs = find(squareGainedInfection == true);

        % Iterate over the newly infected squares.
        for k = 1 : size(newlyInfectedSquareIDs, 2)

            farmSquareID = newlyInfectedSquareIDs(k);

            % Retrieve the row for the farm square from the farm square table.
            farmSquareRow = farmSquareTableMatrix(farmSquareID, :);

            % If the farm square has not been infected previously.
            if farmSquareRow(7) == -1
    
                % Draw the fraction of cattle in the farm square which are dairy cattle
                % from the Uniform distribution.
                dairyCattleFraction = rand * (dairyCattleFractionBounds(2) - ...
                dairyCattleFractionBounds(1)) + dairyCattleFractionBounds(1);
    
                % Instantiate a new Farm_Square_Model object and set it to a demographic
                % steady state with host populations from the farm square table.
                farmSquareModel = Bluetongue_TPT_Model.Initialise_Farm_Square_Model( ...
                farmSquareRow(2), farmSquareRow(3), farmSquareParameters, simulationTime, ...
                round([dairyCattleFraction * farmSquareRow(4), (1 - dairyCattleFraction) * ...
                farmSquareRow(4), farmSquareRow(5)]));
    
                % Update the farm square table to add the time of first infection for the
                % farm square.
                farmSquareTableMatrix(farmSquareID, 7) = simulationTime;
    
                % Add the farm square ID and model object to the infected squares.
                infectedFarmSquareModels{end + 1} = {farmSquareID, farmSquareModel};
                
            % If the farm square has been infected previously.
            else
    
                % Retrieve the Farm_Square_Model object for the farm square
                % from infectedFarmSquareModels.
                for p = 1 : size(infectedFarmSquareModels, 2)
                    if infectedFarmSquareModels{p}{1} == farmSquareID
                        farmSquareModel = infectedFarmSquareModels{p}{2};
                        break;
                    end
                end
    
                fprintf("[INFO] - Re-infection has occurred in farm square %d.\n", farmSquareID);
            end
    
            % Draw the number of infected vectors to insert into the farm square from the
            % Uniform distribution.
            initialInfectedVectorCount = round(rand * (initialInfectedVectorCountBounds(2) - ...
            initialInfectedVectorCountBounds(1)) + initialInfectedVectorCountBounds(1));
    
            % Insert infected vectors to seed infection *within* the farm square.
            farmSquareModel.vectorStates.tensor(:, end, end) = initialInfectedVectorCount;
    
            % Update the farm square table to indicate that the farm square
            % has been (re-)infected, along with the time of infection.
            farmSquareTableMatrix(farmSquareID, 6) = true;
            farmSquareTableMatrix(farmSquareID, 8) = simulationTime;
        end

        fprintf("[INFO] - Between-farm transmission complete. %d new farm square(s) infected." ...
        + " [Time: %0.1f seconds].\n", size(newlyInfectedSquareIDs, 2), ...
        toc(betweenFarmTransmissionStopwatch));

        % ----- END SUBSECTION -----

        %% ----- SUBSECTION: WITHIN-FARM TRANSMISSION -----

        withinFarmTransmissionStopwatch = tic;

        % Construct a vector of the IDs of potentially infected farm squares.
        infectedFarmSquareIDs = zeros(size(infectedFarmSquareModels));

        for i = 1:size(infectedFarmSquareModels, 2)
            infectedFarmSquareIDs(i) = infectedFarmSquareModels{i}{1};
        end

        % Take 'slices' of the sixth and ninth columns of farmSquareTableMatrix so that
        % these can be indexed in the same fashion as infectedFarmSquareIDs in the parfor
        % loop below.
        farmSquareTableMatrixSliceSix = farmSquareTableMatrix(infectedFarmSquareIDs,6);
        farmSquareTableMatrixSliceNine = farmSquareTableMatrix(infectedFarmSquareIDs,9);

        % Initalise a variable to track whether Tau_Leap_Step produces a NaN value in this
        % time step. If it does not, the value remains at 0. If it does, the value is set
        % to the sum of the IDs of the farm squares which threw the error.
        nanErrorThrown = 0;

        % Iterate over the potentially infected farm squares using parallel workers.
        parfor i = 1:size(infectedFarmSquareIDs, 2)

            farmSquareID = infectedFarmSquareIDs(i);

            % Get the Farm_Square_Model object for the potentially infected square.
            farmSquareModel = infectedFarmSquareModels{i}{2};

            % Retrieve the mean daily temperature for the farm square on this day (time
            % value t).
            weatherStationCell = cell2mat(farmSquareStationDictionary(farmSquareID));
            weatherStationID = weatherStationCell(1);
            farmSquareMinMaxTemp = cell2mat(stationTemperatureDictionary(weatherStationID));
            farmSquareTemperature = sum(farmSquareMinMaxTemp(:,simulationTime))/2;
            
            % Perform tau-leap for the within-farm transmission model until we have
            % advanced the within-farm model by tau.
            for j = 1 : 1/tauScale

                exitCode = farmSquareModel.Tau_Leap_Step([1 0 0 0 1 1], tau * tauScale, ...
                farmSquareTemperature);

                if exitCode == 1 && farmSquareTableMatrixSliceSix(i)

                    % In this case, the farm square is supposedly infected but
                    % Tau_Leap_Step has signalled that infection has just ceased
                    % so we update the infection status in the farm square table.
                    farmSquareTableMatrixSliceSix(i) = false;

                elseif exitCode == 2
                    
                    % In this case, all hosts in the farm square have been wiped
                    % out (by disease or otherwise) so we should stop tau-leaping
                    % and update the extinction status in the farm square table.
                    farmSquareTableMatrixSliceNine(i) = true;

                    fprintf("[WARN] Hosts are extinct in farm square %d.\n", farmSquareID);

                    break;

                elseif exitCode == 3
                    nanErrorThrown = nanErrorThrown + farmSquareID;
                    break;
                end
            end

            % We must manually propogate changes to Farm_Square_Model objects back to the
            % client in parfor loops.
            infectedFarmSquareModels{i}{2} = farmSquareModel;
        end

        if nanErrorThrown ~= 0
            input("[ERROR] Tau_Leap_Step has produced at least one NaN value." + ...
            " Something has gone very wrong. Press Return to continue (or use a" + ...
            " breakpoint to pause here). ");
        end

        % Update the farm square table based on changes made during the parfor loop.
        farmSquareTableMatrix(infectedFarmSquareIDs, 6) = farmSquareTableMatrixSliceSix;
        farmSquareTableMatrix(infectedFarmSquareIDs, 9) = farmSquareTableMatrixSliceNine;

        % Add the number of infected farm squares at the end of this time step to
        % infectedFarmTimeSeries.
        infectedFarmTimeSeries(simulationTime - simulationStartTime + 2) = ...
        sum(farmSquareTableMatrix(:,6));

        fprintf("[INFO] - Within-farm transmission (τ-leap) complete for %d farm squares." ...
        + " [Time: %0.1f seconds].\n", size(infectedFarmSquareIDs, 2), ...
        toc(withinFarmTransmissionStopwatch));

        % Increment the simulation time to indicate that a time step has been taken.
        simulationTime = simulationTime + tau;

        % ----- END SUBSECTION -----

        fprintf("[INFO] [Simulation total time: %0.1f seconds].\n", toc(simulationStopwatch));
    end

     % ----- END SECTION -----

     fprintf("\n[INFO] Simulation complete! Total time: %0.1f seconds.\n", toc(simulationStopwatch));
     fprintf("[INFO] Cleaning up...\n");

     % Convert the farm square "table" matrix back to a proper table.
     farmSquareTable = array2table(farmSquareTableMatrix);
     farmSquareTable.Properties.VariableNames = ["ID", "Longitude", "Latitude", ...
     "Cattle Population", "Sheep Population", "Currently Infected?", ...
     "First Infection Time", "Latest Infection Time", "Extinction Occurred?"];
     farmSquareTable = convertvars(farmSquareTable, [6,9], 'logical');

     clear betweenFarmTransmissionStopwatch dairyCattleFraction farmSquareID ...
     farmSquareInfected farmSquareMidMortalityTimeSeries farmSquareModel farmSquareRow ...
     farmSquareTableMatrix farmSquareTableMatrixSliceSix farmSquareTableMatrixSliceNine ...
     generatedValidParams i infectedFarmSquareCache infectedFarmSquareID infectedFarmSquareIDs ...
     initialInfectedVectorCount j k memoryStruct newlyInfectedSquareIDs simulationStopwatch ...
     simulationTime sizeOfInfectedFarmSquareModels squareGainedInfection tau tau_j taudash_j ...
     withinFarmTransmissionStopwatch p farmSquareDistanceMatrix ans;

     %input("[INPUT] Press Return to export the simulation data to a MAT file (or use" ...
     %+ " a breakpoint to pause here). ");

     fprintf("[INFO] Saving data...\n");

     % Export all workspace variables to a MAT file.
     save("-v7.3", sprintf("Outputs/%s.mat", exportFileName));

     % Export first infection times for each farm square to a GeoJSON file.
     Bluetongue_TPT_Model.Export_Infection_GeoJSON("Datasets/GLW3_CATTLE_DA_FRANCE_POINTS.geojson", ...
     farmSquareTable, sprintf("Outputs/%s.geojson", exportFileName));

     % Stop saving Command Window output to a log file.
     diary off;
end