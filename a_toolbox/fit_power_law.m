function [R2, p_values] = fit_power_law(A, B)
    % 幂律模型拟合: B = a * (A + c)^b
    % 输入:
    %   A - 神经数据 (向量)
    %   B - 行为数据 (向量)
    % 输出:
    %   R2 - 拟合优度 (R^2)
    %   p_values - [p_a, p_b, p_c] 参数 a, b, c 的 p 值
    
    % 确保 A 和 B 是列向量
    A = A(:);
    B = B(:);
    
    % 定义幂律函数模型
    powerLaw = @(params, A) params(1) * (A + params(3)).^params(2);
    
    % 初始参数猜测 [a, b, c]
    initialParams = [1, 1, 1];

    % 进行非线性最小二乘拟合
    [paramsFit, resid, J, covB] = nlinfit(A, B, powerLaw, initialParams);

    % 获取拟合参数
    a_fit = paramsFit(1);
    b_fit = paramsFit(2);
    c_fit = paramsFit(3);

    % 计算参数标准误差 (se)
    se = sqrt(diag(covB));

    % 计算 P 值
    t_vals = paramsFit ./ se;  % t 统计量
    df = length(A) - length(paramsFit);  % 自由度
    p_values = 2 * (1 - tcdf(abs(t_vals), df));  % 双侧 t 检验

    % 计算拟合优度 R^2
    B_pred = powerLaw(paramsFit, A);
    SS_res = sum((B - B_pred).^2);  % 残差平方和
    SS_tot = sum((B - mean(B)).^2); % 总平方和
    R2 = 1 - SS_res / SS_tot;  % R^2 计算
end
