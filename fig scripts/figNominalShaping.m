% Show how nominal input shaping results change with different durations

close;
clear;

G = globalData();
m = G.m;
M = G.M;
l = G.l;
g = G.g;

% Input shaping parameters
xEnd = 0.25;    % Desired distance
st = 0.001;     % Simulation time step
omega = sqrt(g / l);
zeta = 0;
Td = 2 * pi / omega;
t1 = 0;                                                            
t2 = Td / 2;
t2 = round(t2, mydigits(st));                                               % Round t2 to nearest time step


% Plotting settings
titleFontSize = 26;
labelFontSize = 24;
axisFontSize = 16;
legendFontSize = 20;

% Color scheme
map = 'turbo';
cmap = colormap(map);
rowmap = 1:round(256/21):256;

% n = 20;
% colormap(turbo(n));
% colorcustom = turbo(n);
% colors = interp1(linspace(-1,1,n), colorcustom, )

% Initialize figure
velocityFig = figure();

% Create linear system: 1D pendulum below cart, no stiffness or damping
A=[0   0  1 0;                                                              % *[x theta x_dot theta_dot]'
   0   0  0 1;
   0   0  0 0;
   0 -g/l 0 0];

% Single input "u" is desired acceleration
B=[   0;                                                                    % *x''_des'
      0;
      1;
    -1/l];

% 4 outputs: cart pos, ball angle, cart vel, & ball angular vel
C = eye(4);

D = zeros(4, 1);

Arigid =   [0 0 1 0;
            0 0 0 1;
            0 0 0 0;
            0 0 0 0];
        
Brigid =   [0;
            0;
            1;
            0];

states = {'x', 'theta', 'x_{dot}', 'theta_{dot}'};
inputs = {'u'};
outputs = {'x', 'theta', 'x_{dot}', 'theta_{dot}'};
sys = ss(A, B, C, D, 'statename', states, 'inputname', inputs, 'outputname', outputs);

% Shape input and simulate motion for different durations
row = 23;    % for array of colors
for tdes = 0.75 * Td : 0.05 : 1.5 * Td
    row = row - 1;
    tf = tdes - t2; % Duration of pre-convolved profile
    t = 0:st:tf;    % Time vector used for convolution
    tc = 0:st:tdes; % Time vector post-convolution
    lentc = length(tc);
    a = xEnd * (60 * (t / tf) - 180 * (t / tf).^2 + 120 * (t / tf).^3) / (tf ^ 2);

    % Generate impulse amplitudes
    K = exp(-zeta * pi / sqrt(1 - zeta ^ 2));     
    A1 = 1 / (1 + K);
    A2 = K / (1 + K);
    
    % Create vector of impulses
    pulses = zeros(1, length(t));                           
    pulses(int16(1 + t1 / st)) = A1;                                        % Impulse 1 at t1
    pulses(int16(1 + t2 / st)) = A2;                                        % Impulse 2 at t2
    pulsetimes = find(pulses);                                              % Get indices of pulses
    pulses = pulses(1:pulsetimes(end));                                     % Trim trailing zeros from vector of pulses

    % Convolve desired trajectory with impulses
    ac = conv(a, pulses);

    % Make time and kinematics arrays same length if rounding error
    if length(ac) > length(tc) % Trim vectors after desired time  
        ac = ac(1:length(tc));
    else % Pad end of array with zeros
        ac(length(tc)) = 0;
    end

    % Simulate using shaped input
    output  = lsim(sys, ac, tc);
 
    % Assign outputs to variables
    vel_sim = output(:, 3);

    % Normalize time
    tc = tc/tc(end);

    % Plot
    set(groot, 'CurrentFigure', velocityFig)
    hold on;
    plot(tc, vel_sim, 'Color', cmap(rowmap(row), :), 'LineWidth', 2)

end

% Post-process plot
set(groot, 'CurrentFigure', velocityFig)
set(gca,'FontSize',axisFontSize)
xlabel('Normalized Time', 'FontSize', labelFontSize)
ylabel('Cup Velocity (m/s)', 'FontSize', labelFontSize)
ylim([0 0.7])

clim([0.75 * Td 1.5 * Td])
colormap turbo
c = colorbar;
c.Label.String = 'Movement Duration (s)';
c.Label.FontSize = labelFontSize;
c.Direction = 'reverse';
c.TickLabels = {'\fontname{Arial}1.45\itT', '\fontname{Arial}1.41\itT', ...
    '\fontname{Arial}1.34\itT', '\fontname{Arial}1.27\itT', ...
    '\fontname{Arial}1.20\itT', '\fontname{Arial}1.13\itT', ...
'\fontname{Arial}1.06\itT', '\fontname{Arial}1.00\itT', ...
'\fontname{Arial}0.92\itT', '\fontname{Arial}0.85\itT', ...
'\fontname{Arial}0.78\itT'};    
c.TickLabelInterpreter = 'tex';

% '0.78T', '0.85T', '0.92T', 'T', ...
%     '1.06T', '1.13T', '1.2T', '1.27T', '1.34T', '1.41T', }