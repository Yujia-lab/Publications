function [r_pos, r_neg, r_com] = a_predict_behavior2(all_vects, all_behav, thresh, fs_option, cov)
% CPM: with leave one out cross-validation
% all_vects: E by K, E being the number of edges, K being the number of subjects
% all_behav: K by 1
% threshold: p-value threshold for edge selection
% fs_option: feature selection options
% 1: Pearson correlation
% 2: Spearman correlation (rank)
% 3: Robust regression
% 4: Partial correlation with a covariate
% cov: covariance K by L when fs_option = 4 to run partial correlation

% all_vects = PSD;
% all_behav = rt_mean_stea';
% thresh = 0.01;

if fs_option == 4
    if nargin < 5
        disp('Need to input the covariate');
        return;
    end
end

[E, no_sub] = size(all_vects);

behav_pred_pos = zeros(no_sub, 1);
behav_pred_neg = zeros(no_sub, 1);
behav_pred = zeros(no_sub, 1);

if fs_option == 1
    disp('Edge selected based on Pearson correlation');
elseif fs_option == 2
    disp('Edge selected based on Spearman correlation');
elseif fs_option == 3
    disp('Edge selected based on robust regression');
elseif fs_option == 4
    disp('Edge selected based on partial correlation');
end

fprintf('Leave one out cross validation');
for leftout = 1:no_sub
    train_vects = all_vects;
    train_vects(:, leftout) = [];
    train_behav = all_behav;
    train_behav(leftout) = [];

    if fs_option == 1
        % correlate all edges with behavior using Pearson correlation
        [r_vec, p_vec] = corr(train_vects', train_behav);
    elseif fs_option == 2
        % correlate all edges with behavior using rank correlation
        [r_vec, p_vec] = corr(train_vects', train_behav, 'type', 'Spearman');
    elseif fs_option == 3
        % correlate all edges with behavior using robust regression
        warning('off');
        r_vec = zeros(E, 1);
        p_vec = zeros(E, 1);
        for edge_i = 1:E
            [~, stats] = robustfit(train_vects(edge_i, :)', train_behav);
            cur_t = stats.t(2);
            r_vec(edge_i) = sign(cur_t) * sqrt(cur_t^2 / (no_sub - 2 + cur_t^2));
            p_vec(edge_i) = 2 * (1 - tcdf(abs(cur_t), no_sub - 2)); % two-tailed
        end
    elseif fs_option == 4
        % correlate all edges with behavior using partial correlation
        cov_train = cov;
        cov_train(leftout, :) = [];
        [r_vec, p_vec] = partialcorr(train_vects', train_behav, cov_train);
    end

    % set threshold and define masks
    pos_mask = (r_vec > 0 & p_vec < thresh);
    neg_mask = (r_vec < 0 & p_vec < thresh);

    %     %-----------------sigmoidal weighting---------------------------%
    % pos_edges = find(r_vec > 0 );
    % neg_edges = find(r_vec < 0 );
    % 
    % % covert p threshold to r threshold
    % T = tinv(thresh/2, no_sub-1-2);
    % R = sqrt(T^2/(no_sub-1-2+T^2));
    % 
    % % create a weighted mask using sigmoidal function
    % % weight = 0.5, when correlation = R/3;
    % % weight = 0.88, when correlation = R;
    % pos_mask(pos_edges) = sigmf( r_vec(pos_edges), [3/R, R/3]);
    % neg_mask(neg_edges) = sigmf( r_vec(neg_edges), [-3/R, R/3]);
    % %---------------sigmoidal weighting-----------------------------%

    % get sum of all edges in TRAIN subs
    train_sumpos = sum(train_vects(pos_mask, :), 1)';
    train_sumneg = sum(train_vects(neg_mask, :), 1)';

    % build model on TRAIN subs
    b = regress(train_behav, [train_sumpos, train_sumneg, ones(no_sub - 1, 1)]);
    b_pos = regress(train_behav, [train_sumpos, ones(no_sub - 1, 1)]);
    b_neg = regress(train_behav, [train_sumneg, ones(no_sub - 1, 1)]);

    % run model on TEST sub
    test_sumpos = sum(all_vects(pos_mask, leftout));
    test_sumneg = sum(all_vects(neg_mask, leftout));

    behav_pred(leftout) = b(1) * test_sumpos + b(2) * test_sumneg + b(3);
    behav_pred_pos(leftout) = b_pos(1) * test_sumpos + b_pos(2);
    behav_pred_neg(leftout) = b_neg(1) * test_sumneg + b_neg(2);
end

[r_pos, p_pos] = corr(behav_pred_pos, all_behav);
[r_neg, p_neg] = corr(behav_pred_neg, all_behav);
[r_com, p_com] = corr(behav_pred, all_behav);

cof_pos = regress(all_behav, [behav_pred_pos,  ones(no_sub, 1)]);
cof_neg = regress(all_behav, [behav_pred_neg, ones(no_sub, 1)]);



error_pos = mean((behav_pred_pos - all_behav).^2);
error_neg = mean((behav_pred_neg - all_behav).^2);
end
