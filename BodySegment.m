classdef BodySegment < handle
    %BodySegment Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        global_proximal_trans   % origin of the segment (proximal is on the pelvis side)
        global_distal_trans     % end-point of the segment
        global_cm_trans         % 
        local_cm_trans          % local transformation from proximal to cm
        local_distal_trans      % local transformation from proximal to distal joint
    end
    
    methods
        function this = BodySegment()
            %BodySegment Construct an instance of this class
            %   Detailed explanation goes here          
        end
        
        function outputArg = method1(this,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = this.Property1 + inputArg;
        end
    end
end

