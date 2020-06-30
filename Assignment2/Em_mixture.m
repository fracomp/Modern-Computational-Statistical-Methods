%The following function begins with a preprocessing step. Then inside the
%loop the M step will perform on the Q function until a stopping criteria
%is met

function [iter,pkp1,sigmakp1] = Em_mixture(file_path,pk,sigmak,maxiter)
    %PRE-PROCESSING STEPS
    %load the file into a variable of str type
    xnp1_str=load(file_path);
    %Transform the str type into a manageable matrix called xnp1
    xnp1=cell2mat(struct2cell(xnp1_str));
    %Create the xn matrix
    xn=circshift(xnp1,1);
    %delete the first observation from both xnp1 and xn
    xnp1=xnp1(2:length(xnp1));
    xn=xn(2:length(xn));
    
    %EM ALGORITHM 
    for iter = 0:maxiter
        
        diff1=(xnp1-(1/2)*xn);
        diff2=(xnp1-xn);
        f1 = (1/sqrt(2*pi))*(exp(-(diff1.^2)./(2.*sigmak.^2))./sigmak);
        f2 = (1/sqrt(2*pi))*(exp(-(diff2.^2)./(2.*sigmak.^2))./sigmak);
        %create/update eik
        ei1=(pk*f1)./(pk*f1 +(1-pk)*f2);
        ei2=((1-pk)*f2)./(pk*f1 +(1-pk)*f2);
        %updating the new value of p and sigma
        pkp1=sum(ei1)/length(xnp1);
        sigmakp1=sqrt(sum(ei1.*(diff1.^2) + ei2.*(diff2.^2))./length(xn));
        
        if all(abs([pkp1-pk; sigmakp1-sigmak])<1e-4)
            break; 
        else
            pk = pkp1; sigmak = sigmakp1;
        end
    end      
end