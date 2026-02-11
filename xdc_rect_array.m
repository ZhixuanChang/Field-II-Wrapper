function aper = xdc_rect_array(elem_pos, elem_width, elem_height, math_size, focus)
% XDC_RECT_ARRAY - 生成任意位置的矩形阵元阵列
%
% 用法:
%   aper = xdc_rect_array(elem_pos, elem_width, elem_height, math_size, focus)
%
% 输入参数:
%   elem_pos    - 阵元中心坐标 [N x 3] 矩阵，每行为 [x, y, z]，单位: m
%   elem_width  - 阵元宽度（x方向），单位: m
%   elem_height - 阵元高度（y方向），单位: m
%   math_size   - 数学子元素尺寸（正方形），单位: m
%   focus       - 固定聚焦点 [x, y, z]，单位: m
%
% 输出参数:
%   aper        - Field II 换能器句柄
%
% 示例:
% 创建8x8平面阵列
%   N_x = 8;
%   N_y = 8;
%   [X, Y] = meshgrid(0:N_x-1, 0:N_y-1);  % 0到7，共8个元素
%   X = (X - (N_x-1)/2) * 0.4e-3;  % 居中并缩放
%   Y = (Y - (N_y-1)/2) * 0.4e-3;
%   elem_pos = [X(:), Y(:), zeros(N_x*N_y, 1)];
%   aper = xdc_rect_array(elem_pos, 0.35e-3, 0.35e-3, 0.1e-3, [0 0 40e-3]);
%
% 作者: AI Assistant
% 日期: 2025

    %% 输入验证
    if nargin < 5
        error('需要5个输入参数');
    end
    
    % 检查 elem_pos 维度
    if size(elem_pos, 2) ~= 3
        error('elem_pos 必须是 N x 3 矩阵');
    end
    
    N_elements = size(elem_pos, 1);  % 阵元数量
    
    % 检查标量参数
    if ~isscalar(elem_width) || ~isscalar(elem_height) || ~isscalar(math_size)
        error('elem_width, elem_height 和 math_size 必须是标量');
    end
    
    % 检查 focus 维度
    if length(focus) ~= 3
        error('focus 必须是 [x, y, z] 向量');
    end
    
    %% 计算每个阵元的数学子元素划分
    % 计算x和y方向的子元素数量
    n_sub_x = ceil(elem_width / math_size);
    n_sub_y = ceil(elem_height / math_size);
    
    % 调整 math_size 以精确填充阵元
    actual_math_width = elem_width / n_sub_x;
    actual_math_height = elem_height / n_sub_y;
    
    %% 生成矩形信息矩阵
    % 每个阵元有 n_sub_x * n_sub_y 个数学子元素
    n_math_per_elem = n_sub_x * n_sub_y;
    total_rects = N_elements * n_math_per_elem;
    
    % 预分配矩阵 (每个矩形19列)
    rect = zeros(total_rects, 19);
    
    rect_idx = 1;  % 当前矩形索引
    
    %% 为每个物理阵元生成数学子元素
    for elem_no = 1:N_elements
        % 当前阵元中心
        elem_center = elem_pos(elem_no, :);
        
        % 生成该阵元的所有数学子元素
        for iy = 1:n_sub_y
            for ix = 1:n_sub_x
                % 计算子元素相对于阵元中心的偏移
                offset_x = (ix - 1 - (n_sub_x - 1) / 2) * actual_math_width;
                offset_y = (iy - 1 - (n_sub_y - 1) / 2) * actual_math_height;
                
                % 子元素中心坐标
                sub_center_x = elem_center(1) + offset_x;
                sub_center_y = elem_center(2) + offset_y;
                sub_center_z = elem_center(3);
                
                % 计算四个角点（顺时针顺序，从左下角开始）
                half_w = actual_math_width / 2;
                half_h = actual_math_height / 2;
                
                % 角点1: 左下 (x-, y-)
                x1 = sub_center_x - half_w;
                y1 = sub_center_y - half_h;
                z1 = sub_center_z;
                
                % 角点2: 右下 (x+, y-)
                x2 = sub_center_x + half_w;
                y2 = sub_center_y - half_h;
                z2 = sub_center_z;
                
                % 角点3: 右上 (x+, y+)
                x3 = sub_center_x + half_w;
                y3 = sub_center_y + half_h;
                z3 = sub_center_z;
                
                % 角点4: 左上 (x-, y+)
                x4 = sub_center_x - half_w;
                y4 = sub_center_y + half_h;
                z4 = sub_center_z;
                
                % 填充矩形信息
                rect(rect_idx, :) = [
                    elem_no,           ... % 物理阵元编号
                    x1, y1, z1,        ... % 角点1
                    x2, y2, z2,        ... % 角点2
                    x3, y3, z3,        ... % 角点3
                    x4, y4, z4,        ... % 角点4
                    1.0,               ... % 加权系数 (默认为1)
                    actual_math_width, ... % 宽度
                    actual_math_height,... % 高度
                    sub_center_x, sub_center_y, sub_center_z  ... % 中心坐标
                ];
                
                rect_idx = rect_idx + 1;
            end
        end
    end
    
    %% 生成物理阵元中心矩阵
    center = elem_pos;
    
    %% 调用 Field II 的 xdc_rectangles 函数
    aper = xdc_rectangles(rect, center, focus);
    
    %% 显示信息
    fprintf('成功创建任意阵列:\n');
    fprintf('  物理阵元数: %d\n', N_elements);
    fprintf('  数学子元素数: %d\n', total_rects);
    fprintf('  每阵元子元素: %d x %d = %d\n', n_sub_x, n_sub_y, n_math_per_elem);
    fprintf('  阵元尺寸: %.3f mm x %.3f mm\n', elem_width*1000, elem_height*1000);
    fprintf('  子元素尺寸: %.3f mm x %.3f mm\n', actual_math_width*1000, actual_math_height*1000);
end