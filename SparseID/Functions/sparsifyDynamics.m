function Xi = sparsifyDynamics(Theta,dXdt,lambda,n)
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz
% 
% if n >1
%     warning('if you are trying to calculate Xi from Theta*Xi = 0 you can only do one variable at a time')
% end
% % compute Sparse regression: sequential least squares
% if dXdt == 0
%     Xinull = null(Theta);
%     ind = find(sum(abs(Theta*Xinull), 1) == min(sum(abs(Theta*Xinull),1)));
%     Xi = Xinull(:, ind);
% else
%     Xi = Theta\dXdt;  % initial guess: Least-squares
% end

 Xi = Theta\dXdt;  % initial guess: Least-squares
 
% lambda is our sparsification knob.
for k=1:10
    smallinds = (abs(Xi)<lambda);   % find small coefficients
    Xi(smallinds)=0;                % and threshold
    for ind = 1:n                   % n is state dimension
        biginds = ~smallinds(:,ind);
        if dXdt ==0
            Xi
            Xi(biginds, ind) = null(Theta(:,biginds));
        else
            % Regress dynamics onto remaining terms to find sparse Xi
            Xi(biginds,ind) = Theta(:,biginds)\dXdt(:,ind)
        end
    end
end