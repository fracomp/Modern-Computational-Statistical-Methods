%Defining a function that accept a filename path, a column where we want to
%read and the initial vector and return the betas estimates

function [beta,k,eucl] = Netwon_Raphson(file_path,column,init)
    %A-Reading the observed values     
    file_all=xlsread(file_path);
    Yi=file_all(:,column);
    
    %B-Creating the X matrix
    x1 = ones(141,1);
    x2 = [1:141]';
    x3 = [cosd(30*x2)];
    x4 = [sind(30*x2)];
    X=[x1,x2,x3,x4];
    
    %C-transposing the initial guess estimate and assign to beta
    beta=init';
    %eucl will measure the distance between 2 successive betas vectors
    eucl=4; 
    % k will count the number of iteration needed to find the best
    % estimates
    k=0;

    %D-do this operation until the beta vector does not change substantially
    while (eucl > 3e-16)
            %Creating mu and M (diag(mu)
            mu= (exp(X*beta));
            M = diag(mu);
            %gradient
            g= X' * (Yi - mu);
            %Hessian
            H=-(X'*(M*X));
            %update beta values
            betak=beta; %this vector keeps the value of bj's at iteration k
            
            beta = beta - (H\g);
            k=k+1;
            eucl=norm(betak - beta);
    end    
end