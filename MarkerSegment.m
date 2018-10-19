classdef MarkerSegment < handle
    %MarkerSegment: Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        markers_names = {}   % cell array of marker names in the segment
        marker_pos           % 3d array of marker positions
    end
    
    methods
        function this = MarkerSegment(markers_names, marker_pos)
            %MarkerSegment Construct an instance of this class
            %   Detailed explanation goes here
            this.markers_names = markers_names;
            this.marker_pos = marker_pos;
        end
        
    end
end

