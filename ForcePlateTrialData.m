classdef ForcePlateTrialData
    %FPlateTrialData Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fplate_trial_file_name   % csv file name of the trial
        fplate_names = {}        % cell array of force plate names
        fplate_units = {}         % cell array of ForcePlateData objs
    end
    
    methods
        function this = ForcePlateTrialData(fplate_trial_file_name)
            %FPlateTrialData Construct an instance of this class
            %   Detailed explanation goes here
            this.fplate_trial_file_name = fplate_trial_file_name;
        end        
    end
end

