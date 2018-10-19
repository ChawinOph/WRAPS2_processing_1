% This function groups the data by time point.

function timedata = sortTime(data, n_markers)

%trim the frame and subframe columns
%data = data(:, 3:3*length(markernames) +2);
[rows columns] = size(data);
timedata = zeros(3, n_markers, rows);

%Sort by frame
for n = 0: n_markers - 1
    timedata(:, n+1, :) = transpose(data(:, 3*n + 1 : 3*n + 3)); 
end

end
