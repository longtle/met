function [T, RC, avgNEigen,pointRC, U_Mat, S_Mat, Error] =  IE_DeltaLam_k_MET(A, k, epsilon)

%   MET algorithm to find k edges whose deletion creates largest drop in leading lambda
%   inputs 
%       - matrix A
%       - budget k
%       - slack variable epsilon
%       
%   outputs
%       - list of edges are deleted: T
%       - staticstics 
%            (i) RC is how many time we recompute
%            (ii)avgNEigen: Average # of eigenvalues tracked. 
%            (iii) the rest are for testing/debuging purpose



%Find the list of edges that we want to remove
T = [];
RC = 0; %Take statics how many times we recompute eigenscores
pointRC = []; % Record when we recompute eigenscores
Error = 0;

nEigen = 2; % Begin with two eigenvalues
trackNEigen = [nEigen];
%Tracking eigenvalues and eigenvectors
U_Mat = [];
S_Mat = [];

previousGap = 1000000; %At beginning, this gap is updated in subsequent recompuations
while (size(T, 1) < k)
    %fprintf ('Eigen-scores round %i \n', RC)
    pointRC = [pointRC, size(T,1)];
    [a,b,c] = find(triu(A,1));
    nEdges = nnz(A)/2;
    
    score = zeros(nEdges, nEigen);
    [U, S, V] = svds(A, nEigen + 1);
    %U_Mat = [U_Mat diag(S)];
    %S_Mat = [S_Mat U];
    RC = RC + 1;
    trackNEigen = [trackNEigen nEigen];
    eigenTrack = diag(S);
    pivot = eigenTrack(end);
    eigenTrack = eigenTrack(1:nEigen);
    
    %Compute eigen-score
    for n=1:nEigen
        %tmpU = abs(U(:, n));
        tmpU = U(:,n);
        score(:, n) = 2* tmpU(a).*tmpU(b);
    end
    
    
    %increase or decrease nEigen
    currentGap = eigenTrack(1) - S(end, end);
    if (currentGap > previousGap)
        newEigen = nEigen - 1;
    end
    if (currentGap < previousGap)
        newEigen = nEigen + 1;
    end
    if (newEigen < 1)
        newEigen = 1;
    end
    previousGap = currentGap;
    
    %Add edge one by one until the estimated largest is higher than
    %pivot value
    [maxEigen, maxEigIdx] = max(eigenTrack);
    E = [];
    while ((maxEigen > (pivot - epsilon)) && (size(T, 1) < k) )
        [~, idx] = max(score(:, maxEigIdx));
        for j = 1:nEigen
            pickScore = score(idx, j);
            eigenTrack(j) = eigenTrack(j) - pickScore;
            score(idx, j) = 0;
        end
        %eigenTrack;
        A(a(idx), b(idx)) = 0;
        A(b(idx), a(idx)) = 0;
        E = [a(idx), b(idx)];
        T = [T; E];
        [maxEigen, maxEigIdx] = max(eigenTrack);
    end
    
    nEigen = newEigen;
    
end
%fprintf ('Number of time we need to re-compute %i \n', RC);
%trackNEigen
%avgNEigen = mean(trackNEigen)-1;

avgNEigen = mean(trackNEigen);
%fprintf ('Track average # eigen %i \n', mean(trackNEigen));
T = T(1:k,:);
