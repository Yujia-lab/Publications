
    function p_value = calculatePValue(SSE, R2, DFE, N)
    % 确定自变量数量
    k = N - DFE - 1; % 自变量数量

    % 计算均方误差
    MSE = SSE / DFE;

    % 计算F统计量
    F = (R2 / k) / ((1 - R2) / DFE);

    % 计算P值
    p_value = 1 - fcdf(F, k, DFE);
end
