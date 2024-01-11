function folder = simulationFolder(count)
       
    % 1. Original Input Shaping (No Impedance External Model)
    if count == 1
        folder = '01. Original Input Shaping'
        
    % 2. Multi-Mode Input Shaping with FF:
    elseif count == 2
        folder = '02. Multi-Mode Internal Model With Feedforward'
        
    % 3. Slow Mode IS with FF:
    elseif count == 3
        folder = '03. Slow Mode Internal Model With Feedforward'
        
    % 4. Fast Mode IS with FF:
    elseif count == 4
        folder = '04. Fast Mode Internal Model With Feedforward'
    
    % 5. Rigid Body IS with FF:
    elseif count == 5
        folder = '05. Rigid Body Internal Model With Feedforward'
    
    % 6. No-Impedance IS with FF:
    elseif count == 6
        folder = '06. No-Impedance Internal Model With Feedforward'   

    % 99. Experimental Trajectories
    elseif count == 99
        % For experimental trajectories folder name doesn't matter. This is
        % used by minVelPlot.m implementation
        folder = 'experiment'
        
    end
    
    folder = strcat(folder, '/');
    
end