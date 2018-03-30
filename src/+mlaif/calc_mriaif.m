% Generate fake MRI curves to be replaced later with real ones
function [mriaif] = calc_mriaif(tn,dt,t1,t2,a1,a2,b1,b2,q2,q3);

import mlaif.*;

% number of time points
N = max(size(tn));

% MRI AIF main component
q1 = 1.0;
n1 = gamma(a1+1.0)/b1^(a1+1.0);
temp1 = q1*dt*(tn.^a1).*exp(-b1*tn)/n1;
for i = 1:N
    if tn(i)<t1
        ca1(i) = 0.0;
        i1 = i;
    else
        ca1(i) = temp1(i-i1);
    end
end

% recirculation part
n2 = gamma(a2+1.0)/b2^(a2+1.0);
temp2 = q2*dt*(tn.^a2).*exp(-b2*tn)/n2;
for i = 1:N
    if tn(i)<t2
        ca2(i) = 0.0;
        i2 = i;
    else
        ca2(i) = temp2(i-i2);
    end
end

% steady state piece
q3 = 0.0015;
sumca = zeros(1,N);
sumca(1) = ca1(1);
for i = 2:N
    sumca(i) = sumca(i-1) + q3*ca1(i);
end

% put all the pieces together
mriaif = ca1 + ca2 + sumca;

return