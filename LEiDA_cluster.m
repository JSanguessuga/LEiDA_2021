function LEiDA_cluster

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LEADING EIGENVECTOR DYNAMICS ANALYSIS
%
% This function computes k-means clustering of leading eigenvectors
%
% - Reads the pre-processed LEiDA_data
% - Saves the Cluster properties for a range of K
% 
% Saves the outputs to LEiDA_k_results.mat
%
% Joana Cabral May 2016
% Modified November 2020
% joana.cabral@med.uminho.pt
% Modified by João Sanguessuga, April 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% USER ADAPT HERE:
load_file = 'LEiDA_EigenVectors'; % file created with LEiDA_data.m
save_file  = 'LEiDA_Kmeans_results';

% Set maximum/minimum number of clusters
mink=2;
maxk=20;

%%

load(load_file,'V1_all');
rangeK = mink:maxk;
Kmeans_results=cell(size(rangeK));

disp('Clustering eigenvectors into:')
    
for K=1:length(rangeK)
    disp(['- ' num2str(rangeK(K)) ' FC states'])
    [IDX, C, SUMD, D]=kmeans(V1_all,rangeK(K),'Distance','Cosine','Replicates',1000,'MaxIter',1000,'Display','final','Options',statset('UseParallel',1));
    [~, ind_sort]=sort(hist(IDX,1:rangeK(K)),'descend');
    [~,idx_sort]=sort(ind_sort,'ascend');
    Kmeans_results{K}.IDX=idx_sort(IDX);     % Cluster time course - numeric collumn vector
    Kmeans_results{K}.C=C(ind_sort,:);       % Cluster centroids (FC patterns)
    Kmeans_results{K}.SUMD=SUMD(ind_sort);   % Within-cluster sums of point-to-centroid distances
    Kmeans_results{K}.D=D(:,ind_sort);       % Distance from each point to every centroid
end

save(save_file,'Kmeans_results','rangeK')

disp(['K-means clustering completed and results saved as ' save_file])

