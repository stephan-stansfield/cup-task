function folder = simulationFolder(count)
       
    % 1. Original Input Shaping (No Impedance Model)
    if count == 1
        folder = "1. Original Input Shaping"
    
    % 2. Multi-Mode Input Shaping, FF:
    elseif count == 2
        folder = "2. Multi-Mode, Feedforward"
        
    % 3. No Impedance IS, FF:
    elseif count == 3
        folder = "3. No Impedance Simplification, Feedforward"
        
    % 4. Rigid Body IS, FF:
    elseif count == 4
        folder = "4. Rigid Body Simplification, Feedforward"
        
    % 5. Slow Mode IS, FF:
    elseif count == 5
        folder = "5. Slow Mode Simplification, Feedforward"
    
    % 6. Fast Mode IS, FF:
    elseif count == 6
        folder = "6. Fast Mode Simplification, Feedforward"
        
    % 7. Multi-Mode Input Shaping, No FF:
    elseif count == 7
        folder = "7. Multi-Mode, No Feedforward"
        
    % 8. No Impedance IS, No FF:
    elseif count == 8
        folder = "8. No Impedance Simplification, No Feedforward"
        
    % 9. Rigid Body IS, No FF:
    elseif count == 9
        folder = "9. Rigid Body Simplification, No Feedforward"
    
    % 10. Slow Mode IS, No FF:
    elseif count == 10
        folder = "10. Slow Mode Simplification, No Feedforward"
    
    % 11. Fast Mode IS, No FF:
    elseif count == 11
        folder = "11. Fast Mode Simplification, No Feedforward"    
        
    % 12. Submovements, No Impedance
    elseif count == 12
        folder = "12. Submovements, No Impedance"
        
    % 13. Submovements, Impedance, Feedforward
    elseif count == 13
        folder = "13. Submovements, Impedance, Feedforward"
        
    % 14. Submovements, Impedance, No Feedforward
    elseif count == 14
        folder = "14. Submovements, Impedance, No Feedforward"
        
    % 15. Rigid Body Simplification, Feedforward, Fixed Impedance            % K = 100, B = 10
    elseif count == 15
        folder = "15. Rigid Body Simplification, Feedforward, Fixed Impedance"
        
    % 16. Rigid Body Simplification, Feedforward, Fixed Impedance, 10k eval, 500ms delay            % K = 100, B = 10
    elseif count == 16
        folder = "16. Rigid Body Simplification, Feedforward, Fixed Impedance, 10k eval, 500ms delay"
        
    % 17. Rigid Body Simplification, Feedforward, B=10, K=75-125, 10k eval, 500ms delay
    elseif count == 17
        folder = "17. Rigid Body Simplification, Feedforward, B=10, K=75-125, 10k eval, 500ms delay"
        
    % 18. Rigid Body Simplification, Feedforward, B=6-10, K=75-125, 10k eval, 500ms delay
    elseif count == 18
        folder = "18. Rigid Body Simplification, Feedforward, B=6-10, K=75-125, 10k eval, 500ms delay"
        
    % 19. Rigid Body Simplification, Feedforward, B=6-10, K=75-125, 10k eval, 500ms delay
    % This simulation calculated duration for the vel-duration correlation
    % based on the simulated duration, not the duration of the experimental
    % trial being fit
    elseif count == 19
        folder = "19. Rigid Body Simplification, Feedforward, B=6-10, K=75-125, 10k eval, 500ms delay"
        
    % 20. Rigid Body Simplification, Feedfwd, B=6-10, K=75-175, 10k eval, +500ms delay
    elseif count == 20
        folder = "20. Rigid Body Simplification, Feedfwd, B=6-10, K=75-175, 10k eval, +500ms delay"
        
    elseif count == 21
        folder = "21. Rigid Body Simplification, Feedfwd, B=6-10, K=75-175, 10k eval, 250-625ms delay"
        
    elseif count == 22
        folder = "22. Rigid Body Simplification, FBO, B=6-10, K=75-175, 10k eval, +500ms delay"
        
    elseif count == 23
        folder = "23. Rigid Body Simplification, Feedforward, updated sim duration"
        
    elseif count == 24
        folder = "24. Rigid Body Simplification, Feedforward, B=6-20, K=0-250, +500ms delay"
        
    elseif count == 25
        folder = "25. Rigid Body Simplification, Feedforward, B=7-20, K=0-250, +500ms delay"
        
    elseif count == 26
%         folder = "26. Original Input Shaping 500ms"
        folder = "1. Original Input Shaping"
        
    elseif count == 27
        folder = "27. Multi-Mode Input Shaping 500ms"
        
    elseif count == 28
        folder = "28. Rigid Body Simplification, Feedforward, B=7-20, K=0-250, +250ms delay"
        
    elseif count == 29
%         folder = "29. Rigid Body Simplification, Feedforward, B=0-50, K=0-750, +500ms delay"
        folder = "28. Rigid Body Simplification, Feedforward, B=7-20, K=0-250, +250ms delay"
        
    elseif count == 30
        folder = "30. Multi-Mode Input Shaping, B=7-20, K=0-250, +250ms"
        
    elseif count == 31
        folder = "31. Original Input Shaping 250ms"
        
    % 99. Experimental Trajectories
    elseif count == 99
        % For experimental trajectories folder name doesn't matter. This is
        % used by minVelPlot.m implementation
        folder = "experiment"
        
    end
    
    folder = strcat(folder, "/");
    
end

%     %%%% Feature-Based Objective Function Results: %%%%
%     % 2a. Multi-Mode Input Shaping, FF:
%     elseif count == 14
%         folder = "2a. Multi-Mode, Feedforward, Feature-Based Objective"
%     
%     % 5a. Slow Mode IS, FF:
%     elseif count == 15
%         folder = "5a. Slow Mode Simplification, Feedforward, Feature-Based Objective"
%     
%     % 6a. Fast Mode IS, FF, Feature-Based Objective:
%     elseif count == 16
%         folder = "6a. Fast Mode Simplification, Feedforward, Feature-Based Objective"
%     
%     %%%% Feature-Based objective Function 2nd Attempt (later points debugged) %%%%
%     % 2b. Multi-Mode Input Shaping, FF:
%     elseif count == 17
%         folder = "2b. Multi-Mode, Feedforward, Feature-Based Objective"
%     
%     % 5b. Slow Mode IS, FF:
%     elseif count == 18
%         folder = "5b. Slow Mode Simplification, Feedforward, Feature-Based Objective"
%     
%     % 6b. Fast Mode IS, FF, Feature-Based Objective:
%     elseif count == 19
%         folder = "6b. Fast Mode Simplification, Feedforward, Feature-Based Objective"
%         
%     %%%% Feature-Based objective Function 3rd Attempt (first 3 changes + ball angle + duration) %%%%
%     % 2c. Multi-Mode Input Shaping, FF:
%     elseif count == 20
%         folder = "2c. Multi-Mode, Feedforward, Feature-Based Objective"
%     
%     % 5c. Slow Mode IS, FF:
%     elseif count == 21
%         folder = "5c. Slow Mode Simplification, Feedforward, Feature-Based Objective"
%     
%     % 6c. Fast Mode IS, FF, Feature-Based Objective:
%     elseif count == 22
%         folder = "6c. Fast Mode Simplification, Feedforward, Feature-Based Objective"
%         
%         
%     %%%% Feature-Based objective Function 4th Attempt (only inter-peak min) %%%%
%     % 2d. Multi-Mode Input Shaping, FF:
%     elseif count == 23
%         folder = "2d. Multi-Mode, Feedforward, Feature-Based Objective"
%     
%     % 5d. Slow Mode IS, FF:
%     elseif count == 24
%         folder = "5d. Slow Mode Simplification, Feedforward, Feature-Based Objective"
%     
%     % 6d. Fast Mode IS, FF, Feature-Based Objective:
%     elseif count == 25
%         folder = "6d. Fast Mode Simplification, Feedforward, Feature-Based Objective"
%         
%     %%%% Feature-Based objective Function 5th Attempt (inter-peak min + residual ball angle) %%%%
%     % 2e. Multi-Mode Input Shaping, FF:
%     elseif count == 26
%         folder = "2e. Multi-Mode, Feedforward, Feature-Based Objective"
%     
%     % 5e. Slow Mode IS, FF:
%     elseif count == 27
%         folder = "5e. Slow Mode Simplification, Feedforward, Feature-Based Objective"
%     
%     % 6e. Fast Mode IS, FF, Feature-Based Objective:
%     elseif count == 28
%         folder = "6e. Fast Mode Simplification, Feedforward, Feature-Based Objective"
%         
%     %%%% Feature-Based objective Function 6th Attempt (first 3 changes, nothing else - 1 subject) %%%%
%     % 2f. Multi-Mode Input Shaping, FF:
%     elseif count == 29
%         folder = "2f. Multi-Mode, Feedforward, Feature-Based Objective"
%     
%     % 5f. Slow Mode IS, FF:
%     elseif count == 30
%         folder = "5f. Slow Mode Simplification, Feedforward, Feature-Based Objective"
%     
%     % 6f. Fast Mode IS, FF, Feature-Based Objective:
%     elseif count == 31
%         folder = "6f. Fast Mode Simplification, Feedforward, Feature-Based Objective"
%         
%     %%%% Feature-Based objective Function 7th Attempt (first 3 changes, nothing else - all subjects) %%%%
%     % 2g. Multi-Mode Input Shaping, FF:
%     elseif count == 32
%         folder = "2g. Multi-Mode, Feedforward, Feature-Based Objective"
%     
%     % 5g. Slow Mode IS, FF:
%     elseif count == 33
%         folder = "5g. Slow Mode Simplification, Feedforward, Feature-Based Objective"
%     
%     % 6g. Fast Mode IS, FF, Feature-Based Objective:
%     elseif count == 34
%         folder = "6g. Fast Mode Simplification, Feedforward, Feature-Based Objective"
%         
%     %%%% Feature-Based objective Function 8th Attempt (all changes, manual edit) %%%%
%     % 2h. Multi-Mode Input Shaping, FF:
%     elseif count == 35
%         folder = "2h. Multi-Mode, Feedforward, Feature-Based Objective"
%     
%     % 5h. Slow Mode IS, FF:
%     elseif count == 36
%         folder = "5h. Slow Mode Simplification, Feedforward, Feature-Based Objective"
%     
%     % 6h. Fast Mode IS, FF, Feature-Based Objective:
%     elseif count == 37
%         folder = "6h. Fast Mode Simplification, Feedforward, Feature-Based Objective"
        