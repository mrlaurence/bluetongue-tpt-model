function Parse_Temperature_Data()
    
    importOptions = detectImportOptions("Datasets\Temperature_Raw\stations_formatted.txt");
    importOptions = setvartype(importOptions, {'Var4','Var5'}, 'char');
    stationTable = readtable("Datasets\Temperature_Raw\stations_formatted.txt", importOptions);
    stationTable = removevars(stationTable, {'Var3'});
    stationTable.Properties.VariableNames = ["Station ID", "Station Name", "Latitude", ...
    "Longitude", "Elevation"];
    
    Latitude = dms2degrees(str2double(split(string(table2array(stationTable(:,3))), ":")));
    Longitude = dms2degrees(str2double(split(string(table2array(stationTable(:,4))), ":")));
    
    stationTable = removevars(stationTable, {'Latitude'});
    stationTable = removevars(stationTable, {'Longitude'});
    
    stationTable = [stationTable(:,1:2), table(Latitude, Longitude), stationTable(:,3:end)];
    
    fileList = dir("Datasets\Temperature_Raw\max_formatted\*");

    tempTable = table;

    for file = fileList(3:end)'

        oneTempTable = readtable(strcat(file.folder, '\', file.name));
    
        idx1 = find(table2array(oneTempTable(:,2)) == 20060101);
        idx2 = find(table2array(oneTempTable(:,2)) == 20091231);
    
        if ~isempty(find(table2array(oneTempTable(idx1:idx2,4) ~= 0)))
            fprintf("WARNING: Temperature data quality check failed. There is/are" ...
            + " %d suspect/missing value(s) in file %s. Parsing anyway...\n", max(size(find(table2array(oneTempTable(idx1:idx2,4) ~= 0)))), file.name);
        end

        if file.name == fileList(3).name
            newTempTable = oneTempTable(idx1:idx2, 2:3);
            newTempTable(:,2) = table(table2array(newTempTable(:,2))/10);
            newTempTable.Properties.VariableNames = ["Date", strcat("Station ", file.name)];
        else
            newTempTable = oneTempTable(idx1:idx2, 3);
            newTempTable(:,1) = table(table2array(newTempTable(:,1))/10);
            newTempTable.Properties.VariableNames = [strcat("Station ", file.name)];
        end

        tempTable = [tempTable, newTempTable];
    end
    save("Station_Maximum_Temperature_Table.mat", "tempTable");

    i = 1;
    while true

        if i > size(stationTable, 1)
            break;
        end

        id = num2str(table2array(stationTable(i, 1)));

        len = strlength(id);
        for j = 1 : 6 - len
            id = strcat('0', id);

        end

        if ~ismember(strcat("Station ", id),tempTable.Properties.VariableNames)
            stationTable = [stationTable(1:i-1,:);stationTable(i+1:end,:)];
        else
            i = i + 1;
        end
    end

    save("Station_Table.mat", "stationTable");
end