classdef MarkerTrialData < handle
    %MarkerTrialData Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        marker_trial_file_name       % csv file name of the trial
        marker_segments = {}         % cell array of MarkerSegmenr objs
    end
    
    methods
        function this = MarkerTrialData(marker_trial_file_name)
            %MarkerTrialData Construct an instance of this class
            %   Detailed explanation goes here
            this.marker_trial_file_name = marker_trial_file_name;
        end
        
    end
end

