%author: Long Le
%Given graph, find the edges selected by MET

clear;

graphDir = './sample-graphs/';

k =1000; % the budget

outputDir = './sample-graphs/';

%Read the input file
edges = csvread(strcat(graphDir, 'sample.csv'));
A  = sparse(edges(:,1), edges(:,2), edges(:,3));

%Edges selected by MET
epsilon = 0.5;
[E, RC, avgNEigen, ~, ~, ~, ~] =  IE_DeltaLam_k_MET(A, k, epsilon);

%Compute eigen-drop after finding the list of removed lisf of edges E
[origLambda, newLambda] =  IE_DeltaLam_GivenT_Simple(A, E);
percentDrop = abs(100*(abs(newLambda(1)) - origLambda(1))/origLambda(1));
fprintf ('Budget: %d \n', k);
fprintf ('Percentage drop in the leading eigenvalue: %.2f \n', percentDrop);