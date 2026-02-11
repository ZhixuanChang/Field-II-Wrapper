function [rfdata, t_start] = calc_rf(fc, fs, tx_aper, rx_aper, tx_delay, tx_apod, scat_pos, scat_amp, c)
% CALC_RF - 计算阵列回波射频数据
%
% 用法:
%   [rfdata, t_start] = calc_rf(fc, fs, tx_aper, rx_aper, tx_delay, ...
%                                tx_apod, scat_pos, scat_amp)
%
% 输入参数:
%   fc        - 中心频率 [Hz]
%   fs        - 采样频率 [Hz]
%   tx_aper   - 发射孔径句柄 (Field II)
%   rx_aper   - 接收孔径句柄 (Field II)
%   tx_delay  - 发射延时矩阵 [nu x nv_tx]，单位: s
%               nu = 发射次数, nv_tx = 发射阵元数
%   tx_apod   - 发射变迹矩阵 [nu x nv_tx]，取值范围: [0, 1]
%   scat_pos  - 散射体位置 [ns x 3]，列为 [x, y, z]，单位: m
%   scat_amp  - 散射体散射系数 [ns x 1] 或 [1 x ns]
%   c         - 介质中的声速 [m/s]
%
% 输出参数:
%   rfdata    - 射频数据 [nt x nv x nu]
%               nt = 时间采样点数
%               nv = 接收阵元数
%               nu = 发射次数
%   t_start   - 起始时间 [s]
%
% 示例:
%   fc = 5e6;
%   fs = 100e6;
%   tx_delay = zeros(10, 64);  % 10次发射，64个阵元
%   tx_apod = ones(10, 64);
%   scat_pos = [0 0 30e-3];
%   scat_amp = 1;
%   [rf, t0] = calc_rf(fc, fs, tx_aper, rx_aper, tx_delay, tx_apod, ...
%                      scat_pos, scat_amp);
%   % rf 的大小为 [nt x 64 x 10]
%
% 作者: AI Assistant
% 日期: 2025

    %% 输入验证
    if nargin < 8
        error('需要8个输入参数');
    end
    
    % 检查 tx_delay 和 tx_apod 维度一致
    if ~isequal(size(tx_delay), size(tx_apod))
        error('tx_delay 和 tx_apod 的维度必须相同');
    end
    
    [nu, nv_tx] = size(tx_delay);  % nu: 发射次数, nv_tx: 发射阵元数
    
    % 检查散射体位置
    if size(scat_pos, 2) ~= 3
        error('scat_pos 必须是 [ns x 3] 矩阵');
    end
    
    ns = size(scat_pos, 1);  % 散射体数量
    
    % 确保 scat_amp 是列向量
    scat_amp = scat_amp(:);
    if length(scat_amp) ~= ns
        error('scat_amp 的长度必须等于散射体数量');
    end
    
    %% 设置 Field II 参数
    set_field('fs', fs);
    set_field('c',  c);
    
    %% 生成激励信号和脉冲响应
    dt = 1 / fs;
    t_pulse = 0:dt:2/fc;  % 2个周期
    
    % 加窗正弦信号
    sig = sin(2*pi*fc*t_pulse);
    sig = sig .* hanning(length(sig))';
    
    % 设置激励和脉冲响应
    xdc_excitation(tx_aper, sig);
    xdc_impulse(tx_aper, sig);
    xdc_impulse(rx_aper, sig);
    
    %% 获取接收阵元数量
    % 通过 xdc_get 获取接收孔径信息
    rx_data = xdc_get(rx_aper, 'rect');
    if isempty(rx_data)
        error('无法获取接收孔径信息，请确保孔径已正确创建');
    end
    nv = length(unique(rx_data(1, :)));  % 物理阵元数（第一行是物理阵元编号）
    
    %% 预分配存储空间
    rf_cells = cell(nu, 1);      % 存储每次发射的 RF 数据
    t_starts = zeros(nu, 1);     % 存储每次发射的起始时间
    
    %% 对每次发射进行仿真
    fprintf('开始仿真 %d 次发射...\n', nu);
    
    for i = 1:nu
        % 设置当前发射的延时
        xdc_focus_times(tx_aper, 0, tx_delay(i, :));
        
        % 设置当前发射的变迹
        xdc_apodization(tx_aper, 0, tx_apod(i, :));
        
        % 计算散射回波（使用 calc_scat_multi 获取所有接收阵元的数据）
        [rf_temp, t_start_temp] = calc_scat_multi(tx_aper, rx_aper, ...
                                                   scat_pos, scat_amp);
        
        % rf_temp 的大小为 [nt_i x nv]，其中 nt_i 可能随发射次数变化
        % 存储结果
        rf_cells{i} = rf_temp;
        t_starts(i) = t_start_temp;
        
        % 显示进度
        if mod(i, max(1, floor(nu/10))) == 0 || i == nu
            fprintf('发射进度: %d/%d (%.1f%%)\n', i, nu, 100*i/nu);
        end
    end
    
    %% 对齐所有 RF 数据
    % 找到最小的起始时间
    t_start = min(t_starts);
    
    % 计算每个 RF 数据需要补零的样本数
    delay_samples = round((t_starts - t_start) * fs);
    
    % 找到最大的数据长度
    max_len = 0;
    for i = 1:nu
        current_len = delay_samples(i) + size(rf_cells{i}, 1);
        if current_len > max_len
            max_len = current_len;
        end
    end
    
    % 预分配输出矩阵 [nt x nv x nu]
    rfdata = zeros(max_len, nv, nu);
    
    % 填充对齐后的数据
    for i = 1:nu
        start_idx = delay_samples(i) + 1;
        nt_i = size(rf_cells{i}, 1);
        end_idx = start_idx + nt_i - 1;
        
        % rf_cells{i} 的大小为 [nt_i x nv]
        rfdata(start_idx:end_idx, :, i) = rf_cells{i};
    end
    
    %% 显示摘要信息
    fprintf('\n仿真完成:\n');
    fprintf('  发射次数: %d\n', nu);
    fprintf('  发射阵元数: %d\n', nv_tx);
    fprintf('  接收阵元数: %d\n', nv);
    fprintf('  散射体数量: %d\n', ns);
    fprintf('  中心频率: %.2f MHz\n', fc/1e6);
    fprintf('  采样频率: %.2f MHz\n', fs/1e6);
    fprintf('  RF 数据大小: [%d x %d x %d]\n', size(rfdata, 1), ...
            size(rfdata, 2), size(rfdata, 3));
    fprintf('  起始时间: %.6f us\n', t_start*1e6);
    fprintf('  时间跨度: %.6f us\n', max_len/fs*1e6);
end