classdef Bluetongue_TPT_Model

    % Static class for running the Bluetongue TPT model across a geographical area,
    % running multiple instances of Farm_Square_Model and handling between-farm
    % transmission. Use Bluetongue_TPT_Model.Run() to run a simulation.
    %
    % AUTHOR: Laurence Dhonau.

    methods(Static)

        % @Bluetongue_TPT_Model/Run.m
        Run();

        function Configure_Figure_Presentation()
            
            % Configures the MATLAB environment to make figures much more readable than
            % the defaults.

            set(groot,"defaultLineLineWidth", 2);
            set(groot,"defaultAxesFontSize", 17);
            set(groot,"defaultTextFontSize", 16);
            set(groot, "defaultAxesColorOrder", [0.1333, 0.5451, 0.1333; ...
            0.6350 0.0780 0.1840; 0.0000 0.4470 0.7410; 0.9290 0.6940 0.1250; ...
            0.4940 0.1840 0.5560; 0.3010 0.7450 0.9330]);
            set(groot, "defaultFigurePosition", [680, 300, 730, 550]);
            set(groot,'defaultAxesXGrid','on');
            set(groot,'defaultAxesYGrid','on');
            set(groot,'defaultAxesXMinorGrid','on','defaultAxesXMinorGridMode','manual');
            set(groot,'defaultAxesYMinorGrid','on','defaultAxesYMinorGridMode','manual');
        end

        function farmMatrix = Import_GLW3_JSON(fileName)

            % Imports a GeoJSON file specifying farm squares in the current directory. The
            % GeoJSON object should be a FeatureCollection, where each Feature is a Point
            % with a property 'VALUE' specifying the number of individuals of a certain
            % type (cattle, sheep, etc.) in the farm square.
            %
            % fileName: Name (or path) of the GeoJSON file to import, including extension.
            %
            % Returns a matrix with 3 columns specifying the longitude, latitude and
            % head count, where each row represents a farm square.

            glwJsonData = jsondecode(char(transpose(fread(fopen(fileName),inf))));

            populationVector = vertcat(vertcat(glwJsonData.features.properties).VALUE);
            
            longLatVector = (horzcat(vertcat(glwJsonData.features.geometry).coordinates))';
            
            farmMatrix = [longLatVector round(populationVector)];

            fclose('all');
        end

        function Construct_Farm_Square_Table()

            % Constructs a table encapsulating information about each farm square in the
            % Bluetongue model. The table has columns: ID, Longitude,
            % Latitude, Cattle Population, Sheep Population, Currently Infected?,
            % First Infection Time, Latest Infection Time, Extinction Occurred?. The ID is
            % generated sequentially and populations are imported from the GLW v3 dataset.
            % Each farm square is initialised as uninfected (first infection time = -1).
            %
            % Saves the resulting table to 'Farm_Square_Table.mat'.

            % Import a matrix specifying longitude, latitude and cattle populations for
            % each farm square in France.
            cattleFarmMatrix = Bluetongue_TPT_Model.Import_GLW3_JSON( ...
            "Datasets/GLW3_CATTLE_DA_FRANCE_POINTS.geojson");

            % Do the above for sheep populations.
            sheepFarmMatrix = Bluetongue_TPT_Model.Import_GLW3_JSON( ...
            "Datasets/GLW3_SHEEP_DA_FRANCE_POINTS.geojson");

            % Sort the matrices by latitude and then longitude (in ascending order).
            cattleFarmMatrix = sortrows(cattleFarmMatrix, [2, 1]);
            sheepFarmMatrix = sortrows(sheepFarmMatrix, [2, 1]);

            % Verify that the latitudes and longitudes of farm squares in the cattle and
            % sheep matrices are the same
            if round(cattleFarmMatrix(:,1:2),5) ~= round(sheepFarmMatrix(:,1:2),5)
                input("[ERROR] The latitudes and longitudes of farm squares for the" + ...
                " sheep and cattles datasets do not match. Something has gone very wrong." + ...
                " Press Return to continue (or use a breakpoint to pause here). ");
            end

            % Compute the number of farm squares.
            farmCount = size(cattleFarmMatrix, 1);

            % Construct the farm square table, generating a numeric ID for each farm
            % square, inserting cattle and sheep populations and initialising 'Currently
            % Infected?' to be false and 'First/Latest Infection Time' to be -1.
            farmSquareTable = table([1:farmCount]', cattleFarmMatrix(:,1), cattleFarmMatrix(:,2), ...
            cattleFarmMatrix(:,3), sheepFarmMatrix(:,3), false(farmCount,1), -1 * ones(farmCount, 1), ...
            -1 * ones(farmCount, 1), false(farmCount,1));

            % Set column names for the farm square table.
            farmSquareTable.Properties.VariableNames = ["ID", "Longitude", "Latitude", ...
            "Cattle Population", "Sheep Population", "Currently Infected?", ...
            "First Infection Time", "Latest Infection Time", "Extinction Occurred?"];

            save("Farm_Square_Table.mat", "farmSquareTable");
        end

        function farmSquareID = Lookup_Farm_Square_Lat_Long(farmTable, latitude, longitude)

            % Get the numeric ID of a farm square in the farm square table, given its
            % latitude and longitude.
            %
            % farmTable: A table encapsulating information about each farm in the
            % Bluetongue model (produced by Bluetongue_TPT_Model.Construct_Farm_Table()).
            % latitude: The latitude of the centre of the farm square (in decimal
            % degrees) to 3 decimal places.
            % longitude: The longitude of the centre of the farm square (in decimal
            % degrees) to 3 decimal places.
            %
            % Returns the ID of the farm square (if it exists), which corresponds to the
            % first column of the farm square table.

            % Extract a 2-column matrix specifying the longitude and latitude of each farm
            % square from the farm square table.
            longLatMatrix = table2array(farmTable(:,2:3));

            [squareExists, farmSquareID] = ismember([longitude latitude], round(longLatMatrix, 3), "rows");
            
            if ~squareExists
                fprintf("[ERROR] Invalid farm square lookup attempted. There is no farm" + ...
                " square with those coordinates.\n");
            end
        end

        function Compute_Lat_Long_Distance_Matrix(latLongMatrix)

            % Computes a matrix giving the pairwise haversine distances (in km) between a
            % collection of points on the Earth's surface. This may take some time to run
            % for a large number of points (around 10 minutes for 10,000 points on an Intel
            % Core i5-11600K).
            %
            % latLongMatrix: A matrix with rows representing points and 2 columns giving
            % the latitude and longitude (in decimal degrees) respectively of the point.
            %
            % Saves the resulting matrix to 'Point_Pairwise_Dist.mat'.
            
            % Compute the number of points
            pointCount = size(latLongMatrix, 1);
            
            % Initalise a pairwise distance matrix
            pointDistanceMatrix = zeros(pointCount, pointCount);
    
            % Compute haversine distances for the upper triangle of the pairwise distance
            % matrix
            for i = 1 : pointCount
                fprintf("[INFO] Computing: i = %d\n", i);
                for j = i : pointCount
                    pointDistanceMatrix(i,j) = haversine(latLongMatrix(i,1:2), latLongMatrix(j,1:2));
                end
            end
            
            % Exploit pairwise symmetry to mirror distances to lower triangle
            pointDistanceMatrix = pointDistanceMatrix + pointDistanceMatrix';
            
            save("Point_Pairwise_Dist.mat", "pointDistanceMatrix");
        end

        function Construct_Station_Temperature_Dictionary()

            % Constructs a dictionary encapsulating information about the minimum/maximum
            % daily temperature at each French weather station in the station table over
            % time. The key gives the station ID and the value is a cell array containing
            % a 2 * 1461 matrix, where the first row gives the minimum daily temperature
            % and the second row the maximum daily temperature at the station for each day
            % between 01/01/2006 and 31/12/2009 inclusive. All temperature data comes from
            % the European Climate and Assessment Dataset.
            %
            % Saves the resulting dictionary to 'Station_Temperature_Dictionary.mat'.

            % Load the station maximum and minimum temperature tables into memory.
            maxTempTableStruct = load("Datasets/Station_Maximum_Temperature_Table.mat");
            minTempTableStruct = load("Datasets/Station_Minimum_Temperature_Table.mat");
            maxTempTable = maxTempTableStruct.tempTable;
            minTempTable = minTempTableStruct.tempTable;
            clear maxTempTableStruct minTempTableStruct;

            % Initalise a new dictionary
            temperatureDictionary = dictionary;

            stationNames = string(maxTempTable.Properties.VariableNames(2:end));

            % Iterate over the columns (stations) in the minimum/maximum temperature
            % tables.
            for i = 1 : size(stationNames, 2)

                % Get row vectors containing the minimum/maximum temperatures for the
                % current station.
                minTempVector = table2array(minTempTable(:, i+1))';
                maxTempVector = table2array(maxTempTable(:, i+1))';
                
                % Parse the numeric station ID from the column name.
                stationNameSplit = strsplit(stationNames(i), " ");
                stationID = str2double(stationNameSplit(2));

                % Add the minimum and maximum temperature series to the dictionary under
                % the station ID.
                temperatureDictionary(stationID) = {[minTempVector; maxTempVector]};
            end

            save("Station_Temperature_Dictionary.mat", "temperatureDictionary");
        end

        function Construct_Vector_Mortality_Dictionary

            % Constructs a dictionary encapsulating information about the vector mortality
            % rate in the BTV model for each 0.1°C temperature increment in the range -50
            % °C to 50 °C inclusive. The key gives the temperature (in multiples of 0.1 °C)
            % and the value gives the mortality rate.
            %
            % Saves the resulting dictionary to 'Vector_Mortality_Dictionary.mat'.

            vectorMortalityDictionary = dictionary;
            
            % Initalise the temperature (measured in 0.1 °C).
            temp = -500;

            % We build the dictionary iteratively (and use a temperature "multiplied by
            % 10") to prevent rounding errors in the dictionary keys)
            while temp <= 500
                vectorMortalityDictionary(temp) = 0.009 * exp(0.16 * temp/10);
                temp = temp + 1;
            end

            save("Vector_Mortality_Dictionary.mat", "vectorMortalityDictionary")
        end

        function parameters = Define_Within_Farm_Parameters()

            % Constructs a structure containing within-farm model parameters. Parameters
            % are drawn from prior distributions.
            %
            % Returns the generated structure.
        
            % Initialise a struct to hold parameters
            parameters = struct();
        
            % ----- TRANSMISSION PARAMETERS -----
            
            % Viral replication rate
            parameters.alpha = normrnd(0.020, 0.0020);
            
            % Threshold temperature for viral replication in vectors
            parameters.T_min = normrnd(13.20, 0.20);
            
            % Probabilities of transmission from vector -> host and host -> vector
            parameters.b = betarnd(21.89, 4.17);
            parameters.beta = betarnd(4.78, 203);
            
            % Feeding preference of vectors for sheep relative to cattle
            parameters.sigma = exprnd(0.095);
            
            % ----- INFECTION PARAMETERS -----
            
            % Reciprocal of mean duration of viraemia in cattle
            parameters.r = 1/normrnd(20.53, 0.90);
        
            % Reciprocal of mean duration of viraemia in sheep
            parameters.rbar = 1/normrnd(16.10, 1.01);
            
            % Number of stages of viraemia in cattle
            parameters.K = ceil(normrnd(5.30, 0.66));
        
            % Number of stages of viraemia in sheep
            parameters.Kbar = ceil(gamrnd(12.28, 1.01));
            
            % Number of stages of vector EIP
            parameters.k = ceil(gamrnd(4.38, 2.51));
            
            % Disease-associated mortality rate in cattle
            parameters.d_C = normrnd(0.0012, 0.00055);
        
            % Disease-associated mortality rate in sheep
            parameters.d_S = normrnd(0.0068, 0.0032);
            
            % ----- VECTOR ACTIVITY PARAMETERS -----
        
            parameters.b_11 = normrnd(-1.59, 0.11);
            parameters.b_21 = normrnd(-3.80, 0.29);
            parameters.b_12 = normrnd(-1.46, 0.071);
            parameters.b_22 = normrnd(-0.99, 0.20);

            parameters.activityNormaliser = 1/max(Farm_Square_Model.Compute_Vector_Activity(parameters, 0:0.01:400, 1));
            
            % ----- TPT PARAMETERS -----
            
            % Reciprocal of mean duration of youth in dairy cattle
            parameters.omega = 1/365;
        
            % Reciprocal of mean duration of youth in beef cattle
            parameters.omegatilde = 1/550;
            
            % Number of stages of youth in dairy cattle
            parameters.L = 3;
        
            % Number of stages of youth in beef cattle
            parameters.Ltilde = 5;
            
            % Reciprocal of mean time till first pregnancy after reaching adulthood for
            % dairy/beef cattle
            parameters.psi = 1/30;
            
            % Reciprocal of mean duration of pregnancy in dairy/beef cattle
            parameters.delta = 1/282;
            
            % Number of stages of pregnancy in dairy/beef cattle
            parameters.M = 5;
            
            % Reciprocal of mean duration of postpartum period in dairy cattle
            parameters.epsilon = 1/45;
        
            % Reciprocal of mean duration of postpartum period in beef cattle
            parameters.epsilontilde = 1/80;
            
            % Probability of successful pregnancy in dairy/beef cattle
            parameters.eta = 0.7;
            
            % Probability of transplacental transmission from an infected dam to her calf
            parameters.zeta = 8/115;
        
            % Proportion of beef calves selected for breeding, rather than culling
            parameters.e = 0.17 / parameters.eta;
        
            % ----- POPULATION PARAMETERS -----
        
            % Mean vector-to-host ratio in the farm square
            parameters.mu_V = gamrnd(6.37, 284);
        
            % Shape parameter for vector-to-host ratio in the farm square
            parameters.s_V = gamrnd(3.31, 0.51);
        
            % Birth/non-disease associated mortality rate for sheep in the farm square
            parameters.b_S = 0.003;
        
            % Cull rate for dairy cattle in postpartum period
            parameters.chi = parameters.epsilon * parameters.eta / 2;
            % parameters.chi = parameters.epsilon / (2 / parameters.eta) * (1/log(2));
        
            % Cull rate for beef cattle in postpartum period
            parameters.chitilde = 0.17 * parameters.epsilontilde;
        end

        function farmSquare = Initialise_Farm_Square_Model(longitude, latitude, parameters, startTime, glwPopulations)

            % Instantiates a new Farm_Square_Model object and sets it to a demographic
            % steady state (ignoring stochastic fluctation) such that the dairy
            % cattle/beef cattle/sheep population sizes equal to those passed in
            % glwPopulations.
            %
            % longitude: The longitude of the centre of the farm square (in decimal
            % degrees).
            % latitude: The latitude of the centre of the farm square (in decimal
            % degrees).
            % parameters: A structure specifying within-farm transmission parameters for
            % use in the model. Bluetongue_TPT_Model.Define_Within_Farm_Parameters()
            % generates a structure of the required form.
            % startTime: The value of time (t, in days) to initialise the farm square
            % model at.
            % glwPopulations: A 3-element row vector with elements specifying the dairy
            % cattle/beef cattle/sheep population in the farm square (likely derived from
            % GLW v3).
            %
            % Returns a Farm_Square_Model with totally susceptible dairy cattle/beef
            % cattle/sheep populations at a demographic steady state with populations
            % equal to those passed in glwPopulations.

            % Define a struct specifying the initial populations to use for a 'steady
            % state finder' simulation. We split the populations from GLW evenly across
            % susceptible compartments since this is generally close to the steady state.
            testPopulations = struct("dairy", table(["X"]', {[glwPopulations(1)]}', {true}'), ...
                                     "beef",  table(["X"]', {[glwPopulations(2)]}', {true}'), ...
                                     "sheep",  table(["X"]', {[glwPopulations(3)]}', {false}'));

            % Instantiate a Farm_Square_Model.
            farmSquare = Farm_Square_Model(latitude, longitude, parameters, 0, testPopulations);

            % Run the farm square simulation for a year (this is enough time to ensure
            % convergence to a demographic steady state) with tau = 1.
            farmSquare.Tau_Leap_Step([1 1 1 1 1 1], 1, 30);

            while(farmSquare.t < 365)
                if farmSquare.Tau_Leap_Step([0 0 0 0 0 0], 1, 30) ~= 1
                    input("[ERROR] Tau_Leap_Step exited with a code other than 1 during farm" + ...
                    " square initialisation. Something has gone very wrong. Press Return to" + ...
                    " continue (or use a breakpoint to pause here). ");
                end
            end

            % Set the host and vector state tensors so they have only one entry in time:
            % the demographic steady state reached.
            farmSquare.dairyCattleStates.tensor = farmSquare.dairyCattleStates.tensor(:,:,end);
            farmSquare.beefCattleStates.tensor = farmSquare.beefCattleStates.tensor(:,:,end);
            farmSquare.sheepStates.tensor = farmSquare.sheepStates.tensor(:,:,end);
            farmSquare.vectorStates.tensor = farmSquare.vectorStates.tensor(:,:,end);

            % Set the farm square time to the requested simulation start time.
            farmSquare.t = startTime;
        end

        function data = Load_Dataset(fileName)

            % Loads into memory a .mat file which contains exactly one variable.
            % 
            % fileName: Name (or path) of the .mat file to import, including extension.
            %
            % Returns the value of the stored variable.

            dataStructure = load(fileName);
            variableName = string(fieldnames(dataStructure));
            data = dataStructure.(variableName);
        end

        function farmSquareMidMortalityTimeSeries = Compute_Mid_Mortality_Time_Series( ...
        farmSquareID, startTime, endTime, farmSquareStationDictionary, ...
        stationTemperatureDictionary, vectorMortalityDictionary)

            % Computes a time series for the vector mortality rate on a given farm square
            % based on the temperature at the weather station closest to the farm square.
            %
            % farmSquareID: The ID of the farm square to compute the time series for.
            % tau_j: The start time for the time series.
            % startTime: The first time to include in the time series, as a non-negative
            % integer number of days, where day 0 is 01/01/2006.
            % endTime: The last time to include in the time series, as a non-negative
            % integer number of days, where day 0 is 01/01/2006.
            % farmSquareStationDictionary: A dictionary containing information on which
            % weather station is closest in distance to each farm square. The file
            % 'Datasets/Farm_Square_Nearest_Station_Dictionary.mat' contains a variable
            % of the required form.
            % stationTemperatureDictionary: A dictionary containing temperature time
            % series for each weather station. The file
            % 'Datasets/Station_Temperature_Dictionary.mat' contains a variable of the
            % required form.
            % vectorMortalityDictionary: A dictionary giving the vector mortality rate for
            % each 0.1°C temperature increment in the range -50 °C to 50 °C inclusive. The
            % file 'Datasets/Vector_Mortality_Dictionary.mat' contains a variable of the
            % required form.
            %
            % Returns the computed time series as a vector.

            % Retrieve the ID of the nearest weather station to the farm square.
            weatherStationCell = cell2mat(farmSquareStationDictionary(farmSquareID));
            weatherStationID = weatherStationCell(1);
            
            % Retrieve a matrix with two rows giving time series for the minimum and maximum
            % daily temperature (in the farm square) respectively.
            farmSquareMinMaxTemp = cell2mat(stationTemperatureDictionary(weatherStationID));
    
            % Compute time series for the minimum and maximum daily vector mortality
            % between tau_j and simulationTime (inclusive) in the farm square.
            farmSquareMinMortalityTimeSeries = vectorMortalityDictionary(10 * ...
            farmSquareMinMaxTemp(1,startTime:endTime));
            farmSquareMaxMortalityTimeSeries = vectorMortalityDictionary(10 * ...
                farmSquareMinMaxTemp(2,startTime:endTime));
    
            % Compute a time series for the midpoint of the minimum and maximum daily vector
            % mortalities between tau_j and simulationTime (inclusive).
            farmSquareMidMortalityTimeSeries =  (farmSquareMinMortalityTimeSeries + ...
            farmSquareMaxMortalityTimeSeries) / 2;
        end

        function Export_Infection_GeoJSON(glwJsonFileName, farmTable)

            % Exports a GeoJSON file encapsulating whether each farm square in France has
            % experienced infection prior to a point in time, and, if so, at what value of
            % time (t) infection was introduced. The GeoJSON object is a FeatureCollection,
            % where each Feature is a Point representing a farm square with a property
            % 'VALUE' specifying the value of t in which infection was introduced to the
            % square (or -1 if the square has never been infected).
            %
            % glwJsonFileName: Name (or path) of a GeoJSON file specifying farm squares
            % in the current directory. The GeoJSON object should be a FeatureCollection,
            % where each Feature is a Point with a property 'VALUE' specifying the number
            % of individuals of a certain type (cattle, sheep, etc.) in the farm square.
            %
            % farmTable: A table encapsulating information about each farm in the
            % Bluetongue model at a point in time (produced by Bluetongue_TPT_Model.
            % Construct_Farm_Table()), from which farm infection status is derived.
            %
            % Saves the resulting file to 'BTV_France_Simulation.geojson'.

            % Import the farm square JSON file to a struct.
            glwJsonStruct = jsondecode(char(transpose(fread(fopen(glwJsonFileName),inf))));

            % Import the Features of the farm square struct to a table.
            glwJsonTable = struct2table(glwJsonStruct.features);

            % Iterate over the farm squares in the table.
            for i = 1 : height(glwJsonTable)

                % Retrieve the longitude and latitude of the farm square.
                squareLongLat = table2array(glwJsonTable(i,3)).coordinates;

                % Lookup the farm square ID using the longitude and latitude.
                squareID = Bluetongue_TPT_Model.Lookup_Farm_Square_Lat_Long(farmTable, ...
                round(squareLongLat(2),3),round(squareLongLat(1), 3));

                % Retrieve the time of first infection from the farm square table (this
                % value is -1 if infection has never occured).
                squareInfecTime = table2array(farmTable(squareID, 7));

                % Replace the VALUE property of the farm square in the Features table with
                % the time of first infection.
                glwJsonTable(i, 2) = table(struct("VALUE", squareInfecTime));
            end

            % Having changed the VALUE properties of each farm square, package the Features
            % table back up into a struct and encode the result to JSON.
            glwJsonStruct.features = table2struct(glwJsonTable);
            outputJson = jsonencode(glwJsonStruct);
    
            % Write the JSON object to BTV_France_Simulation.geojson.
            fprintf(fopen("BTV_France_Simulation.geojson", 'w'), "%s", outputJson);  

            fclose('all');
        end
    end
end