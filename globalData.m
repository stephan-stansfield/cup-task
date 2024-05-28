function G = globalData(name, value)
    
    persistent G_
    
    if isempty(G_)
        G_.m = 1.1; % Mass of ball (kg)
        G_.M = 1.9; % Mass of cup (kg)
        G_.l = 0.5; % Length of the massless pendulum rod (m)
        G_.g = 9.8; % Gravity (m/s^2)
        G_.evalCount = 0;
        G_.velChange_exp = [];
        G_.velPeaks_exp = [];
        G_.velChange_sim = [];
        G_.velPeaks_sim = [];
    end
    
    if nargin > 0
        G_.(name) = value;
    end
    
    G = G_;
    
end