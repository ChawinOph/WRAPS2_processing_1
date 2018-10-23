classdef AnthroModel < handle
    %AnthroModel Summary of this class goes here
    %   Detailed explanation goes here
    properties (Constant)
        g = [0 0 -9.81]'; % gravitaionl acceleration w.r.t. inertial frame
    end
    
    properties
        segments               % cell array of segment objects (can be in differnet class)
        segment_names          % name of the segment
        global_transforms      % transform of all segments (4 x 4 x no_segments)
        total_center_of_mass   % center of mass pos w.r.t inertial frame
    end
    
    methods
        function this = AnthroModel()
            %AnthroModel Construct an instance of this class
            %   Detailed explanation goes here
        end
%         
%         function outputArg = (this, inputArg)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = this.Property1 + inputArg;
%         end
    end
end

