function [origLambda, newLambda] =  IE_DeltaLam_GivenT_Simple(A, T)
%Given list of edges T, we want compute the drop in the leading eigenvalue 
[k, ~] = size(T);

origLambda = eigs(A);

for i=1:k
    e = T(i,:);
    A(e(1), e(2)) = 0;    
    A(e(2), e(1)) = 0;    
end
newLambda = eigs(A);
end
