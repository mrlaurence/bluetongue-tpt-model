function formattedStrings = vecsprintf(formatSpec, formatValues)

    % Uses sprintf to produce n formatted strings with a single format specifier
    % (containing a single formatting operator) and an n-dimensional string array of
    % values to insert into the formatting operator.
    %
    % formatSpec: A string in the form of a formatSpec for sprintf.
    % formatValues: A string array containing values to be inserted into the format
    % operator.
    %
    % Returns a string array of the produced formatted strings. 
    %
    % AUTHOR: Laurence Dhonau.

    formattedStrings = strings(size(formatValues));

    for i = 1 : max(size(formattedStrings))
        formattedStrings(i) = sprintf(formatSpec, formatValues(i));
    end
end