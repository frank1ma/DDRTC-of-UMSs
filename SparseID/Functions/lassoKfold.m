
function [MSEcv,estK]=lassoKfold(X,Y,lambda,kfold,random)

n=size(Y,1);

if random==1
% Data partition
c = cvpartition(size(Y,1),'kfold',kfold);

for i=1:kfold
    % Split data Train and Test set
    idx = test(c,i);
    
    Xtr = X(~idx,:);
    Xte= X(idx,:);
    
    Ytr = Y(~idx,:);
    Yte= Y(idx,:);
    
    % Sparse regression, training through training set
    estK=sparsifyDynamics(Xtr,Ytr,lambda,1);
    
    % error of the k-th fold (test set)
    errorK(i)=((Xte*estK)-Yte)'*((Xte*estK)-Yte);
end


elseif random==0
    % Data split
    lenK=round(size(Y,1)/kfold);
    
    for i=1:kfold
        Xte=X((i-1)*lenK+1:(i*lenK),:);
        Xtr=X;
        Xtr((i-1)*lenK+1:(i*lenK),:)=[];
        
        Yte=Y((i-1)*lenK+1:(i*lenK),:);
        Ytr=Y;
        Ytr((i-1)*lenK+1:(i*lenK),:)=[];
        
    % Sparse regression, training through training set
    estK=sparsifyDynamics(Xtr,Ytr,lambda,1);
    
    % error of the k-th fold (test set)
    errorK(i)=((Xte*estK)-Yte)'*((Xte*estK)-Yte);
    
    end
    
end


MSEcv=sum(errorK)/n;




