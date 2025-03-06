function bestD = evaclusters_elbow(data,k_num)
for k = 2:k_num
    opts = statset('MaxIter', 500, 'Display', 'off');
    [IDX1,C1,sumd1,D1] = kmeans(data,k,'Replicates',5,'options',opts);% kmeans matlab
    [yy,ii] = min(D1');      %% assign points to nearest center

    distort = 0;
    distort_across = 0;
    clear clusts;
    for nn=1:k
        I = find(ii==nn);       %% indices of points in cluster nn
        J = find(ii~=nn);       %% indices of points not in cluster nn
        clusts{nn} = I;         %% save into clusts cell array
        if (length(I)>0)
            mu(nn,:) = mean(data(I,:));               %% update mean
            %% Compute within class distortion
            muB = repmat(mu(nn,:),length(I),1);
            distort = distort+sum(sum((data(I,:)-muB).^2));
            %% Compute across class distortion
            muB = repmat(mu(nn,:),length(J),1);
            distort_across = distort_across + sum(sum((data(J,:)-muB).^2));
        end
    end
    %% Set distortion as the ratio between the within
    %% class scatter and the across class scatter
    distort = distort/(distort_across+eps);

        bestD(k)=distort;
        bestC=clusts;
end
figure; plot(bestD);