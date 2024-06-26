function val = objCalc(output,pos,theta,vel,omega,acc,alpha,weights)
% OBJCALC
%
% Calculates the weighted root mean square error between simulated and
% experimental hand and ball positions, velocities, and accelerations.
    
    % Filter out motion profiles that don't meet criteria.
    % Array of large numbers to replace output of solutions that don't
    % meet criteria. Should be passed over by optimization algorithm.
    nullarray = 100000*ones(length(pos),6);
    
    % If cart velocity solution is ~0 velocity:
    vel_sim     = output(:,3);
    if mean(vel_sim) < 0.05      
        output = nullarray;
    end

    % Assign variable names to simulated output profiles for readability
    pos_sim     = output(:,1);
    theta_sim   = output(:,2);
    vel_sim     = output(:,3);
    omega_sim   = output(:,4);
    acc_sim     = output(:,5);
    alpha_sim   = output(:,6);

    % Only calculate RMSE for times where both motions are defined
    lensim = length(pos_sim);
    lenexp = length(pos);

    if lensim ~= lenexp
        smaller = min(lensim,lenexp);

        % Shorten the longer arrays to agree dimensionally
        if lensim > lenexp
            pos_sim_RMSE     = pos_sim(1:smaller);
            theta_sim_RMSE   = theta_sim(1:smaller);
            vel_sim_RMSE     = vel_sim(1:smaller);
            omega_sim_RMSE   = omega_sim(1:smaller);
            acc_sim_RMSE     = acc_sim(1:smaller);
            alpha_sim_RMSE   = alpha_sim(1:smaller);

            pos_RMSE         = pos;
            theta_RMSE       = theta;
            vel_RMSE         = vel;
            omega_RMSE       = omega;
            acc_RMSE         = acc;
            alpha_RMSE       = alpha;
        else
            pos_RMSE         = pos(1:smaller);
            theta_RMSE       = theta(1:smaller);
            vel_RMSE         = vel(1:smaller);
            omega_RMSE       = omega(1:smaller);
            acc_RMSE         = acc(1:smaller);
            alpha_RMSE       = alpha(1:smaller);

            pos_sim_RMSE     = pos_sim;
            theta_sim_RMSE   = theta_sim;
            vel_sim_RMSE     = vel_sim;
            omega_sim_RMSE   = omega_sim;
            acc_sim_RMSE     = acc_sim;
            alpha_sim_RMSE   = alpha_sim;
        end

    % If no time difference, rename original vectors for RMSE code.
    else
        pos_RMSE         = pos;
        theta_RMSE       = theta;
        vel_RMSE         = vel;
        omega_RMSE       = omega;
        acc_RMSE         = acc;
        alpha_RMSE       = alpha;

        pos_sim_RMSE     = pos_sim;
        theta_sim_RMSE   = theta_sim;
        vel_sim_RMSE     = vel_sim;
        omega_sim_RMSE   = omega_sim;
        acc_sim_RMSE     = acc_sim;
        alpha_sim_RMSE   = alpha_sim;

    end

    % Assign variable names to weights for readability
    wPos    = weights(1);
    wTheta  = weights(2);
    wVel    = weights(3);
    wOmega  = weights(4);
    wAcc    = weights(5);
    wAlpha  = weights(6);
    
    % Difference between simulated & actual profiles at each time
    diffPos     = pos_sim_RMSE   - pos_RMSE';
    diffTheta   = theta_sim_RMSE - theta_RMSE';
    diffVel     = vel_sim_RMSE   - vel_RMSE';
    diffOmega   = omega_sim_RMSE - omega_RMSE';
    diffAcc     = acc_sim_RMSE   - acc_RMSE';
    diffAlpha   = alpha_sim_RMSE - alpha_RMSE';

    % Normalize each variable to total range of experimental trajectory
    rangePos     = max(pos_RMSE)   - min(pos_RMSE);
    rangeTheta   = max(theta_RMSE) - min(theta_RMSE);
    rangeVel     = max(vel_RMSE)   - min(vel_RMSE);
    rangeOmega   = max(omega_RMSE) - min(omega_RMSE);
    rangeAcc     = max(acc_RMSE)   - min(acc_RMSE);
    rangeAlpha   = max(alpha_RMSE) - min(alpha_RMSE);

    % Square the differences
    sqPos   = diffPos.^2;
    sqTheta = diffTheta.^2;
    sqVel   = diffVel.^2;
    sqOmega = diffOmega.^2;
    sqAcc   = diffAcc.^2;
    sqAlpha = diffAlpha.^2;

    % Take RMSE of each variable and normalize to variables range
    % in experimental trial
    RMSEpos = sqrt(wPos*sum(sqPos))/rangePos;
    RMSEtheta = sqrt(wTheta*sum(sqTheta))/rangeTheta;
    RMSEvel = sqrt(wVel*sum(sqVel))/rangeVel;
    RMSEomega = sqrt(wOmega*sum(sqOmega))/rangeOmega;
    RMSEacc = sqrt(wAcc*sum(sqAcc))/rangeAcc;
    RMSEalpha = sqrt(wAlpha*sum(sqAlpha))/rangeAlpha;
    RMSEsum = RMSEpos + RMSEtheta + RMSEvel + RMSEomega + RMSEacc + RMSEalpha;

    val = RMSEsum/(sqrt(wPos)+sqrt(wTheta)+sqrt(wVel)+sqrt(wTheta)+...
        sqrt(wAcc)+sqrt(wAlpha));

end