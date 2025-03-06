function [Call,IDXall] = a_cluster(DATA,cluster_num)
%%% Coding by Aoyujia 2020.7.21
%%% See manu for more instructions

% DATA column: Observation Row: feature

tmp2 = std(DATA);
id_max = DynamicBC_extrema(tmp2);
DATAmx = DATA(:, id_max);
fprintf('K-means Clustering: K=%d\n',cluster_num)
[IDX,C,sumd,D] = kmeans(DATAmx', cluster_num, 'distance', 'sqeuclidean', 'Replicates', 5, 'empty', 'drop');
[IDXall, Call, SUMDall, Dall] = kmeans(DATA', cluster_num, 'distance', 'sqeuclidean', 'Replicates', 1, 'Display', 'iter', 'empty', 'drop', 'Start', C);

