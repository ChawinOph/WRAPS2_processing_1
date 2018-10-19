% sortData imports data from a vicon txt file and separates by marker.
% Right now, I can't figure out how to read the headers correctly, so I
% have to input a list of marker names.

function markerdata = sortMarker(data, n_markers)

[rows columns] = size(data);
markerdata = zeros(rows,3, n_markers);

%sort by marker
for n = 0: n_markers - 1
    markerdata(:,:, n+1) = data(:, 3*n + 1 : 3*n + 3);
end

end