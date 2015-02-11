% Bayes PETMR processing

import mlaif.*;

N = 80;
% get PET data from file and spline and display
load('p7861ho1AS.mat');

dt = 1.0;
tn = 1:N;

iAIF = spline(tCor,cAIF,tn);
iWB = spline(t,WB,tn);
iGM = spline(t,CALLg,tn);
iWM = spline(t,CALLw,tn);

% figure;
% plot(t,WB,t,CALLg,t,CALLw);

% plot the real data
figure;
plot(tn,iAIF,tn,iWB,tn,iGM,tn,iWM);

% Generate fake MRI curves to be replaced later with real ones
NPAR = 12;       % number of parameters
lprob_fun = @logprob_exp;  % set up the lprob function

% initial set up of parameters and models
pflg = zeros(NPAR,1); % -1 regular param, 0 fixed param 
pset = zeros(NPAR,1); % set value for fixed params
pmin = zeros(NPAR,1); % minimal value for param
pmax = zeros(NPAR,1); % maximum value for param
pavg = zeros(NPAR,1); % avergae value for param
psig = zeros(NPAR,1); % std value for param

% set up param 1, ncnt, overall scale
pflg(1) = -1.0;
pmin(1) = 1.0e4;
pavg(1) = 1.0e7;
pmax(1) = 3.0e7;
psig(1) = 0.25*(pmax(1)-pmin(1));
pset(1) = pavg(1);

% set up param 2, cbf, 
pflg(2) = -1.0;
pmin(2) = 0.0001;
pavg(2) = 0.01;
pmax(2) = 0.1;
psig(2) = 0.25*(pmax(2)-pmin(2));
pset(2) = pavg(2);

% set up param 3, cbv, 
pflg(3) = -1.0;
pmin(3) = 0.001;
pavg(3) = 0.1;
pmax(3) = 0.5;
psig(3) = 0.25*(pmax(3)-pmin(3));
pset(3) = pavg(3);

% set up param 4, lambda, pet volume of distribution
pflg(4) = -1.0;
pmin(4) = 0.5;
pavg(4) = 0.9;
pmax(4) = 1.0;
psig(4) = 0.25*(pmax(4)-pmin(4));
pset(4) = pavg(4);

% set up param 5, t1 
pflg(5) = -1.0;
pmin(5) = 5.0;
pavg(5) = 20.0;
pmax(5) = 25.0;
psig(5) = 0.25*(pmax(5)-pmin(5));
pset(5) = pavg(5);

% set up param 6, t2
pflg(6) = -1.0;
pmin(6) = 25.0;
pavg(6) = 30.0;
pmax(6) = 50.0;
psig(6) = 0.25*(pmax(6)-pmin(6));
pset(6) = pavg(6);

% set up param 7, a1 
pflg(7) = -1.0;
pmin(7) = 1.0;
pavg(7) = 6.0;
pmax(7) = 15.0;
psig(7) = 0.25*(pmax(7)-pmin(7));
pset(7) = pavg(7);

% set up param 8, a2
pflg(8) = -1.0;
pmin(8) = 1.0;
pavg(8) = 6.0;
pmax(8) = 20.0;
psig(8) = 0.25*(pmax(8)-pmin(8));
pset(8) = pavg(8);

% set up param 9, b1 
pflg(9) = -1.0;
pmin(9) = 0.5;
pavg(9) = 1.1;
pmax(9) = 4.0;
psig(9) = 0.25*(pmax(9)-pmin(9));
pset(9) = pavg(9);

% set up param 10, b2
pflg(10) = -1.0;
pmin(10) = 0.2;
pavg(10) = 1.1;
pmax(10) = 4.0;
psig(10) = 0.25*(pmax(10)-pmin(10));
pset(10) = pavg(10);

% set up param 11, q2 
pflg(11) = -1.0;
pmin(11) = 0.01;
pavg(11) = 0.1;
pmax(11) = 0.4;
psig(11) = 0.25*(pmax(11)-pmin(11));
pset(11) = pavg(11);

% set up param 12, q3
pflg(12) = -1.0;
pmin(12) = 0.0001;
pavg(12) = 0.0015;
pmax(12) = 0.01;
psig(12) = 0.25*(pmax(12)-pmin(12));
pset(12) = pavg(12);


% Call the MCMC for processing
[parmax avpar] = mcmc(iAIF, iWB, NPAR, lprob_fun, pflg, pset, pmin, pmax, pavg, psig);

% get final plots
[mriaif] = calc_mriaif(tn,dt,avpar(5),avpar(6),avpar(7),avpar(8),avpar(9),avpar(10),avpar(11),avpar(12));
[petaif wbtac] = calc_petaif(mriaif,tn,dt,avpar(1),avpar(2),avpar(3),avpar(4));

figure;
plot(tn,mriaif,tn,petaif);

figure;
plot(tn,iAIF,tn,iWB,tn,petaif,tn,wbtac);

return
