% MCMC_EST:
% Markov Chain Monte-Carlo
% Has the machinery to do a simple Bayesian estimation
%
% NPAR: number of parameters
% xdata: array of independent input data
% ydata: array of dependent input data
% lprob_fun: name of function to calc log prob
% 
% P is a struct that contains:
% pflg: -1 regular param, 0 fixed param 
% pset: set value for fixed params
% pmin: minimal value for param
% pmax: maximum value for param
% pavg: avergae value for param
% psig: std value for param
% pnrm: normalization

function [parmax, avpar] = mcmc_est( xdata, ydata, NPAR, lprob_fun, P )

% constants
NPOP = 50;      % number of population
NPOPREP = 5;    % number of population to replace
NBETA = 50;     % number of temperature steps
NANNEAL = 20;   % number of loops per annealing temp
NPROP = 100;    % number of loops in parameter prob phase
LRG =  1.0e20;

N = max(size(xdata));

% initialize the MCMC
ptmp = zeros(NPAR, 1);
sigma = zeros(NPAR, 1);
parmax = zeros(NPAR, 1);
par = zeros(NPAR, NPOP);
for m = 1:NPOP
    for k = 1:NPAR
        sigma(k) = (P.pmax(k) - P.pmin(k))/10.0;  % why 10???
        if (P.pflg(k) >= 0.0)
            par(k,m) = P.pset(k);
        else
            par(k,m) = P.pavg(k) + P.psig(k)*randn(1,1);
            while (par(k,m)<P.pmin(k) || par(k,m)>P.pmax(k))
                par(k,m) = P.pavg(k) + P.psig(k)*randn(1,1);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% annealing /burn in phase
%%%%%%%%%%%%%%%%%%%%%%%%%%

% loop on the temperature
lpvsb = zeros(NBETA, 1);
parvsb = zeros(NPAR, NBETA);
lprob = zeros(NPOP, 1);
avpar = zeros(NPAR, 1);
sdpar = zeros(NPAR, 1);
parn = zeros(NPAR, 1);
initz = zeros(NPAR, 1);

for i=1:NBETA
    fprintf('temp %d\n', i);
    beta = (1.0/(NBETA-1.0))*i;
    parn = initz;
    lpvsb(i) = 0.0;
            
    % population loop
    for m = 1:NPOP
        ptmp = par(:,m);
        lp0 = lprob_fun( xdata, ydata, NPAR, ptmp, P, beta, 1 );
        lprob(m) = lp0;
            
        for j = 1:NANNEAL
            
            for k = 1:NPAR
                if (P.pflg(k) >= 0.0) 
                    continue; 
                else
                    ptmp = par(:,m);
                    nprop = 0;
                    while (nprop < 50)
                        nprop = nprop + 1;
                        dpar = sigma(k)*randn(1,1);
                        if (ptmp(k)+dpar>=P.pmin(k) && ptmp(k)+dpar<=P.pmax(k)) 
                            break; 
                        end
                    end
                    if (nprop < 50)
                        ptmp(k) = ptmp(k) + dpar;
                    else
                        fprintf('Warning: nprop too large, par %d val %f sigma %f\n',k,ptmp(k),sigma(k));
                    end
                    % new lp
                    lp1 = lprob_fun( xdata, ydata, NPAR, ptmp, P, beta, 1 );
                    
                    % metropolis acceptance
                    dele = lp1 - lp0;
                    if (dele > 0.0 || rand(1,1) < exp(dele))
                        lp0 = lp1;
                        lprob(m) = lp1;
                        parn(k) = parn(k) + 1;
                        par(:,m) = ptmp;
                    end
                end % end pflg if
            end % end k loop
        end % end j loop
    end % end m loop
    
    % get stats
    avpar = initz;
    sdpar = initz;
    for m = 1:NPOP
        ptmp = par(:,m);
        for k = 1:NPAR
            avpar(k) = avpar(k) + ptmp(k);
            sdpar(k) = sdpar(k) + ptmp(k)^2;
        end
        lp1 = lprob_fun( xdata, ydata, NPAR, ptmp, P, beta, 0 );
        lpvsb(i) = lpvsb(i) + lp1/NPOP;
    end
    for k = 1:NPAR
        avpar(k) = avpar(k)/NPOP;
        sdpar(k) = (sdpar(k) - NPOP*avpar(k)^2)/(NPOP - 1);
        if (sdpar(k) > 0.0)
            sdpar(k) = sqrt(sdpar(k));
        else
            sdpar(k) = 0.0;
        end
        parvsb(k,i) = avpar(k);
    end
    
    % adjust proposals, ala Bretthorst
    for k = 1:NPAR
        if (P.pflg(k) >= 0.0) 
            continue; 
        end
        
        if (parn(k) < 0.1*NANNEAL*NPOP)
            sigma(k) = 0.2*sigma(k);
        elseif (parn(k) < 0.2*NANNEAL*NPOP)
            if (sigma(k) > 0.05*sdpar(k))
                sigma(k) = sigma(k)/1.2;
            end
        elseif (parn(k) > 0.8*NANNEAL*NPOP)
            if (sigma(k) < sdpar(k))
                sigma(k) = 5.0*sigma(k);
            end
        elseif (parn(k) > 0.6*NANNEAL*NPOP)
            if (sigma(k) < sdpar(k))
                sigma(k) = 2.0*sigma(k);
            end
        elseif (parn(k) > 0.3*NANNEAL*NPOP)
            if (sigma(k) < sdpar(k))
                sigma(k) = 1.1*sigma(k);
            end
        end
        if (sigma(k) > (P.pmax(k) - P.pmin(k)))
            fprintf('sigma of param %d too large\n',k);
            sigma(k) = P.pmax(k) - P.pmin(k);
        end           
    end        
                
    % print stats every so often
    if (mod(i,10) == 0)
        for k = 1:NPAR
            fprintf('anneal step %d par# %d avg %f std %e proposal %e accept rate %e\n',...
                i,k,avpar(k),sdpar(k),sigma(k),parn(k)/(NANNEAL*NPOP));
        end
    end
    
    % section to replace poor members
    for k = 1:NPOPREP
        lpmax = -LRG;
        lpmin = LRG;
        ilpmax = 0;
        ilpmin = 0;
        for m = 1:NPOP
            if (lprob(m) == 0.0) 
                continue; 
            end
            if (lprob(m) < lpmin)
                lpmin = lprob(m);
                ilpmin = m;
            end
            if (lprob(m) > lpmax)
                lpmax = lprob(m);
                ilpmax = m;
            end
        end
        
        % assure that these are not repeated
        lprob(ilpmin) = 0.0;
        lprob(ilpmax) = 0.0;
        
        par(:,ilpmin) = par(:,ilpmax);
    end
    
end % end i loop on NBETA

% plotting on annealing phase
% for k = 1:NPAR
%     figure(k);
%     plot(1:NBETA, parvsb(k,:));
%     title(['Parameter ',num2str(k),' vs temperature']);
% end
figure;
plot(1:NBETA, lpvsb);
title('Log prob vs temperature');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% proposition/sampling phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% population loop
    parhist = zeros(NPAR, NPOP*NPROP/5);    % for histogram of parameters
    lprobqc = zeros(NPOP, NPROP/5);         % qc on the lprob
    sdeverr = zeros(NPOP*NPROP/5, 1);       % follow the standard dev of error
    for m = 1:NPOP
        ptmp = par(:,m);
        parmax = par(:,m);
        lp0 = lprob_fun( xdata, ydata, NPAR, ptmp, P, 1.0, 1 );
        lprob(m) = lp0;
        lpmax = lp0;
        parn = initz;    
        for j = 1:NPROP
            
            for k = 1:NPAR
                if (P.pflg(k) >= 0.0) 
                    continue; 
                else
                    ptmp = par(:,m);
                    nprop = 0;
                    while (nprop < 50)
                        nprop = nprop + 1;
                        dpar = sigma(k)*randn(1,1);
                        if (ptmp(k)+dpar>=P.pmin(k) && ptmp(k)+dpar<=P.pmax(k)) 
                            break; 
                        end
                    end
                    if (nprop < 50)
                        ptmp(k) = ptmp(k) + dpar;
                    else
                        fprintf('Warning: nprop too large, par %d val %f sigma %f\n',k,ptmp(k),sigma(k));
                    end
                    % new lp
                    lp1 = lprob_fun( xdata, ydata, NPAR, ptmp, P, 1.0, 1 );
                    
                    % store best for future use
                    if (lp1 > lpmax)
                        lpmax = lp1;
                        parmax = ptmp;
                    end
                
                    % metropolis acceptance
                    dele = lp1 - lp0;
                    if (dele > 0.0 || rand(1,1) < exp(dele))
                        lp0 = lp1;
                        lprob(m) = lp1;
                        parn(k) = parn(k) + 1;
                        par(:,m) = ptmp;
                    end
                end % end pflg if
            end % end k loop on NPAR
            
            % histogram parameters 
            if (mod(j,5) == 0)
                lprobqc(m, j/5) = lp0;
                for k = 1:NPAR
                    parhist(k,(m-1)*NPROP/5 + j/5) = ptmp(k);
                end
                lp1 = lprob_fun( xdata, ydata, NPAR, ptmp, P, 1.0, -1 );
                sdeverr((m-1)*NPROP/5 + j/5) = sqrt(lp1/(N-2));
            end
        end % end j loop NPROP
    end % end m loop on NPOP

    % print best fit member of population
    for k = 1:NPAR
        fprintf('best par %d val %f \n',k, parmax(k));
    end
    
    % print final stats
%     avpar = initz;
%     sdpar = initz;
%     for m = 1:NPOP
%         ptmp = par(:,m);
%         for k = 1:NPAR
%             avpar(k) = avpar(k) + ptmp(k);
%             sdpar(k) = sdpar(k) + ptmp(k)^2;
%         end
%     end
%     for k = 1:NPAR
%         avpar(k) = avpar(k)/NPOP;
%         sdpar(k) = (sdpar(k) - NPOP*avpar(k)^2)/(NPOP - 1);
%         if (sdpar(k) > 0.0)
%             sdpar(k) = sqrt(sdpar(k));
%         else
%             sdpar(k) = 0.0;
%         end
%     end
% 
%     for k = 1:NPAR
%         fprintf('par %d avg %f std %f\n', k,avpar(k),sdpar(k));
%     end
    
    avpar = mean(par');
    sdpar = std(par');
    for k = 1:NPAR
        fprintf('par %d avg %f std %f\n', k,avpar(k),sdpar(k));
    end
    
%     xspan = linspace(min(ydata(:)), 1.01*max(ydata(:)), NPAR); %NONC
%     figure;
%     title('Final histogram');
%     errorbar(xspan,avpar,sdpar);
    
    % histogram sampling phase
    % histogram parameter distribution
    for k = 1:NPAR
        figure;
        hist(parhist(k, :),50);
        title(['Parameter ',num2str(k),' Distribution']);
    end
    
    % histogram std error
    figure;
    hist(sdeverr, 50);
    title('STD of error Distribution');
    
    % histogram QC of log probability
    figure;
    hold on;
    xprop = 1:NPROP/5; 
    
    for k = 1:NPOP
        plot(xprop, lprobqc(k,:));
    end
    title('Log Prob of Population');
    hold off;
    
end % close function
    






