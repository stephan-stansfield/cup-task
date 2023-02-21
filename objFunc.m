function [val, gradient] = objFunc(x,optimizationType,ampMod,...
        forwardF,intModel,impedance,tc,tdes,delayMin,delayMax,xEnd,vStart,vEnd,...
        st,shift,simVersion,pendIndex,objective,pos,vel,...
        acc,theta,omega,alpha,weights,blockStr,blockNum,subjNum,num,tDesMax,fitMethod)
    
    printDes = false;

    % Increment evaluation counter
%     G = globalData();
%     evalCount = G.evalCount;
%     % DEBUG
%     disp(evalCount)
%     evalCount(num+1) = evalCount(num+1) + 1;
%     G = globalData('evalCount', evalCount);
    

    % Parameters to optimize over are passed in array x. For readability, 
    % assign variables to x array elements.
    b = x(1);
    k = x(2);
    if optimizationType == "input shaping 4 impulse"
        p       = x(3);
        q       = x(4);
        r       = x(5);
        s       = x(6);
        tdelay  = x(7);

        if ampMod == true
            fa11    = x(8);
            fa12    = x(9);
            fa21    = x(10);
            fa22    = x(11);
        else
            fa11    = 1;
            fa12    = 1;
            fa21    = 1;
            fa22    = 1;
        end

    elseif optimizationType == "input shaping 2 impulse impedance"
        p       = x(3);
        q       = x(4);
        tdelay  = x(5);

        if ampMod == true
            fa1   = x(6);
            fa2   = x(7);
        else
            fa1   = 1;
            fa2   = 1;
        end

    elseif optimizationType == "input shaping 2 impulse no impedance"
        p       = x(1);                                                 % In no impedance model, overwrite b and k parameters    
        q       = x(2);
        tdelay  = x(3);

        if ampMod == true
            fa1   = x(4);
            fa2   = x(5);
        else
            fa1   = 1;
            fa2   = 1;
        end

    elseif optimizationType == "submovement"
        D1      = x(3);
        tf1     = x(4);
        ti2     = x(5);
        tdelay  = x(6);
    end

    % Create external system model with chosen parameters.
    % - "print" input is set to "false" because this will get run every 
    % evaluation in optimization
    [extSys,intSys,sysRigid,Td,Td1,Td2,zeta,zeta1,zeta2,overdamped] = ...
        sysCreate(b,k,forwardF,intModel,impedance,false);                

    % Array of large numbers to replace output of solutions that don't
    % meet criteria. Should be passed over by optimization algorithm.
    nullarray = 100000*ones(length(tc),6);

    % Add from 0 up to delayMax seconds to simulated time duration
    % TIME DELAY CHANGE
    tdessim = tdes + delayMin + 0.5*(tdelay-1)*(delayMax-delayMin);         % range from delayMin to delayMax
%     tdessim = tdes + (tdelay-2)*delayMax;                                 % range from -delayMax to +delayMax (MARK FOR DELETION)

    % If system is not overdamped, simulate motion using given inputs
    if optimizationType == "input shaping 4 impulse"

        if overdamped
            output = nullarray;
        else

            modes = 2;
            [output, ~, ~, ~, ~] = simInputShape(b, k, intModel, extSys, ...
                sysRigid, Td1, Td2, zeta1, zeta2, tdes, tdessim, xEnd, ...
                vStart, vEnd, st, shift, forwardF, simVersion, modes, ...
                pendIndex, fitMethod);

        end

    elseif optimizationType == "input shaping 2 impulse impedance"
        
        if overdamped

            output = nullarray;

        else

            modes = 1;
            [output, ~, ~, ~, ~] = simInputShape(b, k, intModel, extSys, ...
                sysRigid, Td, Td2, zeta, zeta2, tdes, tdessim, xEnd, ...
                vStart, vEnd, st, shift, forwardF, simVersion, modes, ...
                pendIndex, fitMethod);
        end

    elseif optimizationType == "input shaping 2 impulse no impedance"

        modes   = 1;
        [output, ~, ~, ~, ~] = simInputShape(b, k, intModel, extSys, ...
            sysRigid, Td, Td2, zeta, zeta2, tdes, tdessim, xEnd, vStart, ...
            vEnd, st, shift, forwardF, simVersion, modes, pendIndex, ...
            fitMethod);

    elseif optimizationType == "submovement"

        [output, ~] = simSubmovements(extSys,sysRigid,tdes,tdessim,st,xEnd,...
        vStart,D1,tf1,ti2,pendIndex,simVersion,forwardF,impedance);

    end

    % Now that motion was simulated, calculate objective function value
    val = objCalc(objective,output,pos,theta,vel,omega,acc,alpha,weights,...
        blockStr,blockNum,subjNum,num);

end