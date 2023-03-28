function minMat = Find_Nearest_Station(stationLatLongMatrix, farmLatLongMatrix)

    farmCount = size(farmLatLongMatrix, 1);
    stationCount = size(stationLatLongMatrix, 1);

    pointDistanceMatrix = zeros(size(farmLatLongMatrix, 1), size(stationLatLongMatrix,1));

    for i = 1 : farmCount
        fprintf("Computing: i = %d\n", i);
        for j = 1 : stationCount
            pointDistanceMatrix(i,j) = haversine(farmLatLongMatrix(i,1:2), stationLatLongMatrix(j,1:2));
        end
    end
    
    minMat = zeros(farmCount, 2);

    for k = 1 : farmCount
        smol = min(pointDistanceMatrix(k,:));
        minMat(k,1) = smol;
        minMat(k,2) = find(pointDistanceMatrix(k,:) == smol);
    end
end