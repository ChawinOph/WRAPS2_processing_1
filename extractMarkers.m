function [var, var_indice] = extractMarkers(sorted_marker_names, sbj_name ,headernames, M)
var_indice = [];
for i = 1:length(sorted_marker_names)
    for j = 1:size(headernames, 2)
        if strcmp(headernames{j}, [sbj_name, ':', sorted_marker_names{i}])
            var_indice = [var_indice, j]; %#ok<AGROW>
        end
    end
end
var = [];
for i = 1:length(var_indice)
    var = [var, M(:, 3*(var_indice(i)-1) + 1 : 3*(var_indice(i)-1) + 3)]; %#ok<AGROW>
end
end

