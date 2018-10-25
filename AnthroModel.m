classdef AnthroModel < handle
    %AnthroModel Summary of this class goes here
    %   Detailed explanation goes here
    properties (Constant)
        g = [0 0 -9.81]'; % gravitaionl acceleration w.r.t. inertial frame
    end
    
    properties        
        global_transforms           % transform of all segments (4 x 4 x no_segments)
        total_center_of_mass        % center of mass pos w.r.t inertial frame
        segment_parameters          % name of the segments, required markers
        segment_transforms          % transformations of each segment in each trial
        raw_marker_data             % transformations of    
    end
    
    methods
        function this = AnthroModel()
            % AnthroModel Construct an instance of this class
            %   Detailed explanation goes here
        end
       
        function
            
            
        end
        

    end
end

