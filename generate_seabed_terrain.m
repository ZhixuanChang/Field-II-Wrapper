function terrain_matrix = generate_seabed_terrain(dx, h, mean_depth)
% GENERATE_SEABED_TERRAIN 随机生成一维局部微观海底地形线
%
% 输入参数:
%   dx - 水平方向上的步进距离 (分辨率，单位：米)
%   h  - 地形的尺度范围/总长度 (单位：米)
%   mean_depth - 平均水深（米）
%
% 输出参数:
%   terrain_matrix - [M, 2] 的矩阵。
%                    其中 M = round(h/dx)
%                    第 1 列为 x 坐标 (水平距离)
%                    第 2 列为 z 坐标 (海底深度)

    % 1. 计算所需的数据点数 M
    M = round(h / dx);

    % 2. 生成 x 坐标向量 (确保为 M 行 1 列的列向量)
    x = linspace(0, h, M)';

    % --- 内部地形特征参数 (可根据仿真需求微调) ---
    amplitude = 0.25;   % 起伏幅度标准差 (米，控制地形粗糙高差)
    beta = 2.2;         % 谱指数 (2.2 适合泥沙或小碎石底质)
    % ---------------------------------------------

    % 3. FFT 长度优化：取大于 M 的下一个 2 的幂次方 (N >= M)
    % 这样不仅 FFT 计算最快，还能有效避免直接使用首尾相接造成的周期性边界假象
    N = 2^nextpow2(M);
    if N < 512
        N = 512; % 保证样本基数足够大，以维持良好的分形统计特性
    end

    % 4. 生成高斯白噪声
    white_noise = randn(1, N);

    % 5. 傅里叶变换到频域
    F = fft(white_noise);

    % 6. 构造频率向量与 1/f^(beta/2) 滤波器
    freq = [0:(N/2-1), -N/2:-1];
    freq(1) = 1e-6; % 将零频(直流)分量设为极小值，防止除以 0
    
    filter_response = (abs(freq)).^(-beta/2);
    F_filtered = F .* filter_response;

    % 7. 逆傅里叶变换回到时域
    topo_raw = real(ifft(F_filtered));

    % 8. 截取我们需要的前 M 个点
    topo_raw = topo_raw(1:M);
    
    % 强制转换为列向量 (M 行 1 列)
    topo_raw = topo_raw(:);

    % 9. 标准化并映射到实际的深度与幅度
    topo_normalized = (topo_raw - mean(topo_raw)) / std(topo_raw); 
    z = topo_normalized * amplitude + mean_depth;         

    % 10. 组合成 [M, 2] 矩阵输出
    terrain_matrix = [x, z];
end