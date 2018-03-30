% function to calculate the PET time courses from MRI and param data
function [petaif wbtac] = calc_petaif(mriaif,tn,dt,ncnt,cbf,cbv,lambda)

    import mlaif.*;

    % number of time points
    N = max(size(tn));

    % PET AIF from MRI
    expMRI2PET = dt*cbf*(1/cbv-1/lambda)*exp(-cbf*(1/cbv-1/lambda)*tn);
    temp1 = conv(mriaif,expMRI2PET);
    petaif = temp1(1:N);

    % PET residue curves
    exp1 = dt*cbf*exp(-(cbf/lambda)*tn);

    % convolution to get WB pet signal
    temp2 = conv(petaif, exp1);
    wbtac = temp2(1:N);

    % display or not AIF and WB activity curves
%     if (1)
%         figure;
%         hold on;
%         plot(tn,mriaif,tn,petaif,'Linewidth',2);
%         title('AIF MRI and PET');
%         hold off;
% 
%         figure;
%         hold on;
%         plot(tn,petaif,tn,wbtac,'Linewidth',2);
%         title('PET AIF and WB tac');
%         hold off;
%     end
    
    petaif = ncnt*petaif;
    wbtac = ncnt*wbtac;
    
    
return
    






