function eva = a_evacluster(DATA,eva_num)

%% Estimate optimal cluster number

tmp2 = std(DATA);
id_max = DynamicBC_extrema(tmp2);
DATAmx = DATA(:, id_max);

IDXest = zeros(size(DATAmx,2),eva_num);
distortion= [];
for i=1:eva_num
    [IDXest(:,i),Cest{i},sumdest{i}] = kmeans(DATAmx',i,'emptyaction','singleton','replicate',5, 'empty', 'drop');
    distortion(i,1) = sum(sumdest{i});
end

variance=distortion(1:end-1)-distortion(2:end);
distortion_percent=cumsum(variance)/(distortion(1)-distortion(end));

% plot
crit = {'silhouette','CalinskiHarabasz', 'DaviesBouldin'};    %     crit = {'CalinskiHarabasz', 'DaviesBouldin', 'gap', 'silhouette'};
fig = figure('color','w','units','norm','pos',[0.3,0.3,0.6,0.3]);
%https://en.wikipedia.org/wiki/Determining_the_number_of_clusters_in_a_data_set
subplot(1,4,1); plot(2:eva_num,distortion_percent,'b-o')
xlabel('Number of Clusters');
ylabel('Percent of variance explained') %Percentage of variance explained is the ratio of the between-group variance to the total variance,
OptimalK=[];

for i=1:length(crit)
    eva{i} = evalclusters(DATAmx',IDXest,crit{i});
    subplot(1,4,i+1); plot(eva{i})%         subplot(2,2,i); plot(eva{i})
    tmp=[];
    tmp.OptimalK = eva{i}.OptimalK;
    OptimalK(i) = eva{i}.OptimalK;
    fprintf('(%s) Optimal K=%d\n',crit{i},OptimalK(i))
    tmp.OptimalY = eva{i}.OptimalY;
    tmp.InspectedK = eva{i}.InspectedK;
    tmp.CriterionValues = eva{i}.CriterionValues;
    tmp.CriterionName = eva{i}.CriterionName;
    tmp.ClusteringFunction = eva{i}.ClusteringFunction;
    tmp.NumObservations = eva{i}.NumObservations;
    tmp.Missing = eva{i}.Missing;
    eva{i} = tmp;
end