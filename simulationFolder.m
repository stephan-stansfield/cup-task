function folder = simulationFolder(count)
% SIMULATION FOLDER
%
% Simple dictionary for assigning names of simulation types.
       
    if count == 1
        folder = '01. Original Input Shaping'
    elseif count == 2
        folder = '02. Multi-Mode Internal Model With Feedforward'
    elseif count == 3
        folder = '03. Slow Mode Internal Model With Feedforward'
    elseif count == 4
        folder = '04. Fast Mode Internal Model With Feedforward'
    elseif count == 5
        folder = '05. Rigid Body Internal Model With Feedforward'
    elseif count == 6
        folder = '06. No-Impedance Internal Model With Feedforward'   
    elseif count == 99
        % For experimental trajectories folder name doesn't matter. This is
        % used by minVelPlot.m implementation
        folder = 'experiment'
    end
    
    folder = strcat(folder, '/');
    
end