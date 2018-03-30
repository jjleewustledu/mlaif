% logprob_exp 
% returns the log probability of y = a*exp(-b*x)
%
% flag: -1  return Q
%       0   return log of normalized Q
%       1   return with beta and priors
%
function lprob = logprob_exp(iAIF,iWB,NPAR, par, pflg, pavg, psig, beta, flag )


    import mlaif.*;

    PARPEN = 0.0; % -1.0 for parameter penalty
    N = max(size(iAIF));
    dt = 1.0;
    tn = 1:N;
    
    % NOTE 
    [mriaif] = calc_mriaif(tn,dt,par(5),par(6),par(7),par(8),par(9),par(10),par(11),par(12));
    [petaif wbtac] = calc_petaif(mriaif,tn,dt,par(1),par(2),par(3),par(4));
    
    derr = (iAIF - petaif).^2 + (iWB - wbtac).^2;
    lprob = sum(derr(:));
    
    if (flag == -1)
        return
    end
    
    % use for sigma fit
    if (par(1) < 0.0)
        lprob = lprob / (-2.0*par(1)^2);
        lprob = lprob + (-0.5*N)*log(2.0*pi*par(1)^2);
    else  % no sigma, use t distribution Jeffrey Prior
        lprob = -0.5*N*log(0.5*lprob);
    end

    % add in beta and pretest probabilities 
    if (flag == 1)
        % beta only operates on likelihoods
        lprob = beta*lprob;
        for i = 1:NPAR
            if (pflg(i) < 0.0)
                lprob = lprob + (PARPEN - (par(i) - pavg(i))^2/(2.0*psig(i)^2));
            end
        end
    end
    
    return
        
end

