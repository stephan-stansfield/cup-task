function [val, gradient] = objFunc(x,optimizationType,forwardF,intModel,...
        impedance,tc,tdes,delayMin,delayMax,xEnd,vStart,st,simVersion,...
        pendIndex,pos,vel,acc,theta,omega,alpha,weights)
% OBJFUNC
%
% Creates a system using a set of hyperparameters from the optimization
% algorithm (by calling sysCreate.m), simulates the system using input
% shaping (by calling simInputShape.m), and returns the objective function
% value (calculated by objCalc.m).

    % Parameters to optimize over are passed in array x. For readability, 
    % assign x array elements to variables.
    if optimizationType == "input shaping 2 impulse no impedance"
        b = 0;
        k = 0;
        tdelay = x(1);
    else
        b = x(1);
        k = x(2);
        tdelay = x(3);
    end

    % Create external system model with chosen parameters. "print" input
    % is set to false because this will get run every evaluation in
    % optimization
    [extSys,~,sysRigid,Td,Td1,Td2,zeta,zeta1,zeta2,overdamped] = ...
        sysCreate(b,k,intModel,impedance,false);                

    % Array of large numbers to replace output of solutions that don't
    % meet criteria. Should be passed over by optimization algorithm.
    nullarray = 100000*ones(length(tc),6);

    % Add from 0 up to delayMax seconds to simulated time duration
    tdessim = tdes + delayMin + 0.5*(tdelay-1)*(delayMax-delayMin);

    % If system is not overdamped, simulate motion using given inputs
    if optimizationType == "input shaping 4 impulse"
        if overdamped
            output = nullarray;
        else
            modes = 2;
            [output, ~, ~] = simInputShape(b, k, intModel, extSys, ...
                sysRigid, Td1, Td2, zeta1, zeta2, tdes, tdessim, xEnd, ...
                vStart, st, forwardF, simVersion, modes, ...
                pendIndex, impedance);
        end
    elseif optimizationType == "input shaping 2 impulse impedance"
        if overdamped
            output = nullarray;
        else
            modes = 1;
            [output, ~, ~] = simInputShape(b, k, intModel, extSys, ...
                sysRigid, Td, Td2, zeta, zeta2, tdes, tdessim, xEnd, ...
                vStart, st, forwardF, simVersion, modes, ...
                pendIndex, impedance);
        end
    elseif optimizationType == "input shaping 2 impulse no impedance"
        modes = 1;
        [output, ~, ~] = simInputShape(b, k, intModel, extSys, ...
            sysRigid, Td, Td2, zeta, zeta2, tdes, tdessim, xEnd, vStart, ...
            st, forwardF, simVersion, modes, pendIndex, impedance);
    end

    % Now that motion was simulated, calculate objective function value
    val = objCalc(output,pos,theta,vel,omega,acc,alpha,weights);
end