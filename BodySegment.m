classdef BodySegment < handle
    %BodySegment Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        global_proximal_trans       % origin of the segment (proximal is on the pelvis side)
        global_distal_trans         % end-point of the segment
        global_cm_trans             % 
        local_cm_trans              % local transformation from proximal to cm
        local_distal_trans          % local transformation from proximal to distal joint
        vicon_segment_names         % used vicon segment names in raw data
        vicon_marker_names          % rows of used vicon marker names. Each row is based on each segment
    end
    
    methods
        function this = BodySegment()
            %BodySegment Construct an instance of this class
            %   Detailed explanation goes here          
        end
       
    end
end

