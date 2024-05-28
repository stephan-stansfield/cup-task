function [posNorm, thetaNorm, velNorm, omegaNorm, accNorm, alphaNorm]...
        = normalize(pos,theta,vel,omega,acc,alpha,tdes,st,stNorm,tc)
% NORMALIZE
%
% Takes 6 kinematic trajectory variables and interpolates them to be the
% same length.
   
    stScale = stNorm*tdes;
    tScale = 0:stScale:tdes;
    
    % Interpolate data so each trial has the same number of points
    if isempty(tc)
        t           = 0:st:tdes;
    else
        t = tc;
    end
    posNorm     = interp1(t,pos,tScale);
    thetaNorm   = interp1(t,theta,tScale);
    velNorm     = interp1(t,vel,tScale);
    omegaNorm   = interp1(t,omega,tScale);
    accNorm     = interp1(t,acc,tScale);
    alphaNorm   = interp1(t,alpha,tScale);
    
    % Recursively replace any NaN value with the value in the index
    % immediately preceding it
    while ~isempty(find(isnan(posNorm),1))
        nanIndices = find(isnan(posNorm));
        posNorm(nanIndices) = posNorm(nanIndices - 1);
    end
    
    while ~isempty(find(isnan(velNorm),1))
        nanIndices = find(isnan(velNorm));
        velNorm(nanIndices) = velNorm(nanIndices - 1);
    end
    
    while ~isempty(find(isnan(accNorm),1))
        nanIndices = find(isnan(accNorm));
        accNorm(nanIndices) = accNorm(nanIndices - 1);
    end
    
    while ~isempty(find(isnan(thetaNorm),1))
        nanIndices = find(isnan(thetaNorm));
        thetaNorm(nanIndices) = thetaNorm(nanIndices - 1);
    end
    
    while ~isempty(find(isnan(omegaNorm),1))
        nanIndices = find(isnan(omegaNorm));
        omegaNorm(nanIndices) = omegaNorm(nanIndices - 1);
    end
    
    while ~isempty(find(isnan(alphaNorm),1))
        nanIndices = find(isnan(alphaNorm));
        alphaNorm(nanIndices) = alphaNorm(nanIndices - 1);
    end
    
end
