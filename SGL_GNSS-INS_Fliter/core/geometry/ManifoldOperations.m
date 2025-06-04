classdef ManifoldOperations < handle
% MANIFOLDOPERATIONS 复合流形操作管理类
%
% 实现复合流形 M = R³×R³×SO(3)×R⁶×R² 上的操作
% 包括流形操作符 ⊞, ⊟, ⊕ 和相关的几何计算
%
% 状态向量结构: X = [pos(3), vel(3), att_log(3), ba(3), bg(3), clk(2)]
% 其中 att_log 是SO(3)在李代数so(3)中的对数表示
%
% 作者: Y
% 日期: 2025.6

    properties (Constant)
        % 状态向量各分量的索引
        POS_IDX = 1:3;      % 位置 R³
        VEL_IDX = 4:6;      % 速度 R³
        ATT_IDX = 7:9;      % 姿态 so(3) (李代数表示)
        BA_IDX = 10:12;     % 加速度计偏差 R³
        BG_IDX = 13:15;     % 陀螺仪偏差 R³
        CLK_IDX = 16:17;    % 时钟参数 R²
        
        STATE_DIM = 17;     % 总状态维度
        TANGENT_DIM = 15;   % 切空间维度 (3+3+3+3+3+0=15, SO(3)在切空间为3维)
    end
    
    methods (Static)
        
        function X_new = boxplus(X, xi)
            % BOXPLUS 复合流形上的加法操作符 ⊞
            %
            % 输入:
            %   X - 17x1 状态向量 (流形上的点)
            %   xi - 15x1 扰动向量 (切空间中的扰动)
            %
            % 输出:
            %   X_new - 17x1 更新后的状态向量
            %
            % 实现文章公式(5): X ⊞ ξ
            
            if length(X) ~= ManifoldOperations.STATE_DIM
                error('状态向量维度错误，期望 %d，实际 %d', ManifoldOperations.STATE_DIM, length(X));
            end
            
            if length(xi) ~= ManifoldOperations.TANGENT_DIM
                error('切向量维度错误，期望 %d，实际 %d', ManifoldOperations.TANGENT_DIM, length(xi));
            end
            
            X_new = zeros(ManifoldOperations.STATE_DIM, 1);
            
            % 位置更新 (R³上的加法)
            X_new(ManifoldOperations.POS_IDX) = X(ManifoldOperations.POS_IDX) + xi(1:3);
            
            % 速度更新 (R³上的加法)
            X_new(ManifoldOperations.VEL_IDX) = X(ManifoldOperations.VEL_IDX) + xi(4:6);
            
            % 姿态更新 (SO(3)上的乘法)
            % 当前姿态的旋转矩阵
            R_current = ManifoldOperations.log_to_rotation(X(ManifoldOperations.ATT_IDX));
            % 应用扰动
            R_new = SO3.boxplus(R_current, xi(7:9));
            % 转换回对数表示
            X_new(ManifoldOperations.ATT_IDX) = ManifoldOperations.rotation_to_log(R_new);
            
            % IMU偏差更新 (R³上的加法)
            X_new(ManifoldOperations.BA_IDX) = X(ManifoldOperations.BA_IDX) + xi(10:12);
            X_new(ManifoldOperations.BG_IDX) = X(ManifoldOperations.BG_IDX) + xi(13:15);
            
            % 时钟参数更新 (R²上的加法) - 没有对应的切向量扰动
            X_new(ManifoldOperations.CLK_IDX) = X(ManifoldOperations.CLK_IDX);
        end
        
        function xi = boxminus(X2, X1)
            % BOXMINUS 复合流形上的减法操作符 ⊟
            %
            % 输入:
            %   X2, X1 - 17x1 状态向量 (流形上的点)
            %
            % 输出:
            %   xi - 15x1 差值向量 (切空间中的差值)
            %
            % 实现文章公式(5): X2 ⊟ X1
            
            if length(X1) ~= ManifoldOperations.STATE_DIM || length(X2) ~= ManifoldOperations.STATE_DIM
                error('状态向量维度错误');
            end
            
            xi = zeros(ManifoldOperations.TANGENT_DIM, 1);
            
            % 位置差值 (R³上的减法)
            xi(1:3) = X2(ManifoldOperations.POS_IDX) - X1(ManifoldOperations.POS_IDX);
            
            % 速度差值 (R³上的减法)
            xi(4:6) = X2(ManifoldOperations.VEL_IDX) - X1(ManifoldOperations.VEL_IDX);
            
            % 姿态差值 (SO(3)上的对数差)
            R1 = ManifoldOperations.log_to_rotation(X1(ManifoldOperations.ATT_IDX));
            R2 = ManifoldOperations.log_to_rotation(X2(ManifoldOperations.ATT_IDX));
            xi(7:9) = SO3.boxminus(R2, R1);
            
            % IMU偏差差值 (R³上的减法)
            xi(10:12) = X2(ManifoldOperations.BA_IDX) - X1(ManifoldOperations.BA_IDX);
            xi(13:15) = X2(ManifoldOperations.BG_IDX) - X1(ManifoldOperations.BG_IDX);
        end
        
        function X_new = oplus(X, v)
            % OPLUS 外部输入操作符 ⊕
            %
            % 对于GNSS/INS系统，外部输入通常作用在切空间上
            % 因此 ⊕ 操作与 ⊞ 操作相同
            %
            % 输入:
            %   X - 17x1 状态向量
            %   v - 外部输入向量 (维度可变)
            %
            % 输出:
            %   X_new - 更新后的状态向量
            
            if length(v) == ManifoldOperations.TANGENT_DIM
                % 切空间扰动
                X_new = ManifoldOperations.boxplus(X, v);
            else
                error('外部输入维度不匹配: %d', length(v));
            end
        end
        
        function J = jacobian_boxplus(X, xi)
            % JACOBIAN_BOXPLUS 计算 ⊞ 操作关于扰动的雅可比矩阵
            %
            % 输入:
            %   X - 17x1 状态向量
            %   xi - 15x1 扰动向量
            %
            % 输出:
            %   J - 17x15 雅可比矩阵 ∂(X⊞ξ)/∂ξ
            %
            % 用于误差状态线性化，实现文章公式(6)
            
            J = zeros(ManifoldOperations.STATE_DIM, ManifoldOperations.TANGENT_DIM);
            
            % 位置部分 (单位矩阵)
            J(ManifoldOperations.POS_IDX, 1:3) = eye(3);
            
            % 速度部分 (单位矩阵)
            J(ManifoldOperations.VEL_IDX, 4:6) = eye(3);
            
            % 姿态部分 (左雅可比矩阵)
            if nargin > 1 && length(xi) >= 9
                xi_att = xi(7:9);
            else
                xi_att = zeros(3, 1);
            end
            
            % 当前姿态
            omega_current = X(ManifoldOperations.ATT_IDX);
            
            % 左雅可比矩阵 (用于SO(3)的切空间映射)
            J_left = SO3.left_jacobian(xi_att);
            
            % 姿态雅可比矩阵的计算比较复杂，这里使用数值微分
            J(ManifoldOperations.ATT_IDX, 7:9) = ManifoldOperations.attitude_jacobian_numerical(X);
            
            % IMU偏差部分 (单位矩阵)
            J(ManifoldOperations.BA_IDX, 10:12) = eye(3);
            J(ManifoldOperations.BG_IDX, 13:15) = eye(3);
            
            % 时钟参数部分没有对应的切向量
        end
        
        function J = jacobian_boxminus(X1, X2)
            % JACOBIAN_BOXMINUS 计算 ⊟ 操作的雅可比矩阵
            %
            % 输入:
            %   X1, X2 - 17x1 状态向量
            %
            % 输出:
            %   J - 结构体，包含对X1和X2的雅可比矩阵
            %     J.wrt_X1 - 15x17 矩阵，∂(X2⊟X1)/∂X1
            %     J.wrt_X2 - 15x17 矩阵，∂(X2⊟X1)/∂X2
            
            J = struct();
            
            % 对X1的雅可比 (负号)
            J.wrt_X1 = zeros(ManifoldOperations.TANGENT_DIM, ManifoldOperations.STATE_DIM);
            J.wrt_X1(1:3, ManifoldOperations.POS_IDX) = -eye(3);
            J.wrt_X1(4:6, ManifoldOperations.VEL_IDX) = -eye(3);
            J.wrt_X1(7:9, ManifoldOperations.ATT_IDX) = -ManifoldOperations.attitude_jacobian_numerical(X1);
            J.wrt_X1(10:12, ManifoldOperations.BA_IDX) = -eye(3);
            J.wrt_X1(13:15, ManifoldOperations.BG_IDX) = -eye(3);
            
            % 对X2的雅可比 (正号)
            J.wrt_X2 = zeros(ManifoldOperations.TANGENT_DIM, ManifoldOperations.STATE_DIM);
            J.wrt_X2(1:3, ManifoldOperations.POS_IDX) = eye(3);
            J.wrt_X2(4:6, ManifoldOperations.VEL_IDX) = eye(3);
            J.wrt_X2(7:9, ManifoldOperations.ATT_IDX) = ManifoldOperations.attitude_jacobian_numerical(X2);
            J.wrt_X2(10:12, ManifoldOperations.BA_IDX) = eye(3);
            J.wrt_X2(13:15, ManifoldOperations.BG_IDX) = eye(3);
        end
        
        function R = log_to_rotation(omega_log)
            % LOG_TO_ROTATION 从李代数表示转换为旋转矩阵
            %
            % 输入:
            %   omega_log - 3x1 李代数向量
            %
            % 输出:
            %   R - 3x3 旋转矩阵
            
            R = SO3.exp(omega_log);
        end
        
        function omega_log = rotation_to_log(R)
            % ROTATION_TO_LOG 从旋转矩阵转换为李代数表示
            %
            % 输入:
            %   R - 3x3 旋转矩阵
            %
            % 输出:
            %   omega_log - 3x1 李代数向量
            
            omega_log = SO3.log(R);
        end
        
        function valid = is_valid_state(X)
            % IS_VALID_STATE 检查状态向量的有效性
            %
            % 输入:
            %   X - 17x1 状态向量
            %
            % 输出:
            %   valid - 逻辑值，状态是否有效
            
            if length(X) ~= ManifoldOperations.STATE_DIM
                valid = false;
                return;
            end
            
            % 检查是否包含NaN或Inf
            if any(~isfinite(X))
                valid = false;
                return;
            end
            
            % 检查姿态部分的有效性
            try
                R = ManifoldOperations.log_to_rotation(X(ManifoldOperations.ATT_IDX));
                valid = SO3.is_valid_rotation(R);
            catch
                valid = false;
            end
        end
        
        function X_corrected = enforce_constraints(X)
            % ENFORCE_CONSTRAINTS 强制满足流形约束
            %
            % 输入:
            %   X - 17x1 状态向量 (可能违反约束)
            %
            % 输出:
            %   X_corrected - 17x1 修正后的状态向量
            
            X_corrected = X;
            
            % 修正姿态部分 - 确保对应有效的旋转矩阵
            try
                R = ManifoldOperations.log_to_rotation(X(ManifoldOperations.ATT_IDX));
                if ~SO3.is_valid_rotation(R)
                    R_corrected = SO3.nearest_rotation(R);
                    X_corrected(ManifoldOperations.ATT_IDX) = ...
                        ManifoldOperations.rotation_to_log(R_corrected);
                end
            catch
                % 如果转换失败，设置为单位旋转
                X_corrected(ManifoldOperations.ATT_IDX) = zeros(3, 1);
            end
            
            % 其他约束检查 (根据物理意义)
            % 例如：限制偏差的合理范围
            max_bias_acc = 10; % m/s²
            max_bias_gyro = 10 * pi / 180; % rad/s
            
            X_corrected(ManifoldOperations.BA_IDX) = ...
                max(-max_bias_acc, min(max_bias_acc, X_corrected(ManifoldOperations.BA_IDX)));
            X_corrected(ManifoldOperations.BG_IDX) = ...
                max(-max_bias_gyro, min(max_bias_gyro, X_corrected(ManifoldOperations.BG_IDX)));
        end
        
        function dist = distance(X1, X2)
            % DISTANCE 计算复合流形上两点间的距离
            %
            % 输入:
            %   X1, X2 - 17x1 状态向量
            %
            % 输出:
            %   dist - 标量，测地距离
            
            % 使用加权距离
            weights = [1, 1, 1, 0.1, 0.1, 0.1, 10, 10, 10, 1, 1, 1, 10, 10, 10]; % 各分量权重
            
            % 计算各分量距离
            pos_dist = norm(X2(ManifoldOperations.POS_IDX) - X1(ManifoldOperations.POS_IDX));
            vel_dist = norm(X2(ManifoldOperations.VEL_IDX) - X1(ManifoldOperations.VEL_IDX));
            
            % 姿态距离 (测地距离)
            R1 = ManifoldOperations.log_to_rotation(X1(ManifoldOperations.ATT_IDX));
            R2 = ManifoldOperations.log_to_rotation(X2(ManifoldOperations.ATT_IDX));
            att_dist = SO3.distance(R1, R2);
            
            ba_dist = norm(X2(ManifoldOperations.BA_IDX) - X1(ManifoldOperations.BA_IDX));
            bg_dist = norm(X2(ManifoldOperations.BG_IDX) - X1(ManifoldOperations.BG_IDX));
            clk_dist = norm(X2(ManifoldOperations.CLK_IDX) - X1(ManifoldOperations.CLK_IDX));
            
            % 加权总距离
            dist = sqrt(weights(1)^2 * pos_dist^2 + weights(4)^2 * vel_dist^2 + ...
                       weights(7)^2 * att_dist^2 + weights(10)^2 * ba_dist^2 + ...
                       weights(13)^2 * bg_dist^2 + weights(16)^2 * clk_dist^2);
        end
        
        function [state_names, tangent_names] = get_state_names()
            % GET_STATE_NAMES 获取状态变量名称
            %
            % 输出:
            %   state_names - 17x1 cell数组，状态变量名称
            %   tangent_names - 15x1 cell数组，切向量变量名称
            
            state_names = {
                'pos_x', 'pos_y', 'pos_z', ...           % 位置
                'vel_x', 'vel_y', 'vel_z', ...           % 速度
                'att_log_x', 'att_log_y', 'att_log_z', ... % 姿态(李代数)
                'ba_x', 'ba_y', 'ba_z', ...              % 加速度计偏差
                'bg_x', 'bg_y', 'bg_z', ...              % 陀螺仪偏差
                'clk_bias', 'clk_drift'                  % 时钟参数
            };
            
            tangent_names = {
                'delta_pos_x', 'delta_pos_y', 'delta_pos_z', ...     % 位置扰动
                'delta_vel_x', 'delta_vel_y', 'delta_vel_z', ...     % 速度扰动
                'delta_att_x', 'delta_att_y', 'delta_att_z', ...     % 姿态扰动
                'delta_ba_x', 'delta_ba_y', 'delta_ba_z', ...        % 加速度计偏差扰动
                'delta_bg_x', 'delta_bg_y', 'delta_bg_z'             % 陀螺仪偏差扰动
            };
        end
        
        function X_init = create_initial_state(pos, vel, att_euler, ba, bg, clk)
            % CREATE_INITIAL_STATE 创建初始状态向量
            %
            % 输入:
            %   pos - 3x1 初始位置 [m]
            %   vel - 3x1 初始速度 [m/s]
            %   att_euler - 3x1 初始姿态欧拉角 [rad] (ZYX顺序)
            %   ba - 3x1 加速度计偏差 [m/s²]
            %   bg - 3x1 陀螺仪偏差 [rad/s]
            %   clk - 2x1 时钟参数 [偏差m, 漂移m/s]
            %
            % 输出:
            %   X_init - 17x1 初始状态向量
            
            % 默认值
            if nargin < 1, pos = zeros(3, 1); end
            if nargin < 2, vel = zeros(3, 1); end
            if nargin < 3, att_euler = zeros(3, 1); end
            if nargin < 4, ba = zeros(3, 1); end
            if nargin < 5, bg = zeros(3, 1); end
            if nargin < 6, clk = zeros(2, 1); end
            
            % 欧拉角转旋转矩阵再转李代数
            R_init = SO3.from_euler(att_euler, 'ZYX');
            att_log = ManifoldOperations.rotation_to_log(R_init);
            
            % 组装状态向量
            X_init = [pos(:); vel(:); att_log(:); ba(:); bg(:); clk(:)];
        end
        
        function [pos, vel, att_euler, ba, bg, clk] = extract_state_components(X)
            % EXTRACT_STATE_COMPONENTS 提取状态向量的各个分量
            %
            % 输入:
            %   X - 17x1 状态向量
            %
            % 输出:
            %   pos - 3x1 位置 [m]
            %   vel - 3x1 速度 [m/s]
            %   att_euler - 3x1 姿态欧拉角 [rad]
            %   ba - 3x1 加速度计偏差 [m/s²]
            %   bg - 3x1 陀螺仪偏差 [rad/s]
            %   clk - 2x1 时钟参数 [m, m/s]
            
            pos = X(ManifoldOperations.POS_IDX);
            vel = X(ManifoldOperations.VEL_IDX);
            
            % 李代数转欧拉角
            R = ManifoldOperations.log_to_rotation(X(ManifoldOperations.ATT_IDX));
            att_euler = SO3.to_euler(R, 'ZYX');
            
            ba = X(ManifoldOperations.BA_IDX);
            bg = X(ManifoldOperations.BG_IDX);
            clk = X(ManifoldOperations.CLK_IDX);
        end
        
        function display_state(X, name)
            % DISPLAY_STATE 显示状态向量内容 (调试用)
            %
            % 输入:
            %   X - 17x1 状态向量
            %   name - 状态名称 (可选)
            
            if nargin < 2
                name = 'State Vector';
            end
            
            [pos, vel, att_euler, ba, bg, clk] = ManifoldOperations.extract_state_components(X);
            
            fprintf('\n=== %s ===\n', name);
            fprintf('位置 [m]: [%8.3f, %8.3f, %8.3f]\n', pos);
            fprintf('速度 [m/s]: [%8.3f, %8.3f, %8.3f]\n', vel);
            fprintf('姿态 [deg]: [%8.2f, %8.2f, %8.2f]\n', rad2deg(att_euler));
            fprintf('加速度计偏差 [m/s²]: [%8.4f, %8.4f, %8.4f]\n', ba);
            fprintf('陀螺仪偏差 [deg/s]: [%8.4f, %8.4f, %8.4f]\n', rad2deg(bg));
            fprintf('时钟参数 [m, m/s]: [%8.3f, %8.6f]\n', clk);
            fprintf('状态有效性: %s\n', logical_to_string(ManifoldOperations.is_valid_state(X)));
            fprintf('=======================\n');
        end
        
    end
    
    methods (Static, Access = private)
        
        function J_att = attitude_jacobian_numerical(X)
            % ATTITUDE_JACOBIAN_NUMERICAL 数值计算姿态部分的雅可比矩阵
            %
            % 对于SO(3)到so(3)的映射雅可比，使用数值微分方法
            %
            % 输入:
            %   X - 17x1 状态向量
            %
            % 输出:
            %   J_att - 3x3 姿态雅可比矩阵
            
            eps = 1e-8;
            J_att = zeros(3, 3);
            
            omega0 = X(ManifoldOperations.ATT_IDX);
            
            for i = 1:3
                % 正向扰动
                omega_plus = omega0;
                omega_plus(i) = omega_plus(i) + eps;
                
                % 负向扰动
                omega_minus = omega0;
                omega_minus(i) = omega_minus(i) - eps;
                
                % 中心差分
                J_att(:, i) = (omega_plus - omega_minus) / (2 * eps);
            end
        end
        
    end
    
    methods (Static, Hidden)
        
        function run_tests()
            % RUN_TESTS 运行复合流形操作的单元测试
            
            fprintf('=== 复合流形操作单元测试 ===\n');
            
            % 测试1: 流形操作符的逆性
            fprintf('测试1: 流形操作符逆性...');
            X1 = ManifoldOperations.create_initial_state(...
                [100; 200; 50], [5; 0; -2], [0.1; 0.2; 0.3], ...
                [0.01; -0.02; 0.005], [0.001; 0.002; -0.001], [10; 0.1]);
            X2 = ManifoldOperations.create_initial_state(...
                [105; 198; 48], [4.8; 0.2; -2.1], [0.15; 0.18; 0.35], ...
                [0.015; -0.018; 0.008], [0.0015; 0.0018; -0.0012], [12; 0.05]);
            
            xi = ManifoldOperations.boxminus(X2, X1);
            X2_recovered = ManifoldOperations.boxplus(X1, xi);
            error1 = norm(X2 - X2_recovered);
            
            if error1 < 1e-10
                fprintf(' 通过 (误差: %.2e)\n', error1);
            else
                fprintf(' 失败 (误差: %.2e)\n', error1);
            end
            
            % 测试2: 状态有效性检查
            fprintf('测试2: 状态有效性检查...');
            valid1 = ManifoldOperations.is_valid_state(X1);
            
            % 创建无效状态
            X_invalid = X1;
            X_invalid(7) = NaN; % 姿态分量设为NaN
            valid2 = ManifoldOperations.is_valid_state(X_invalid);
            
            if valid1 && ~valid2
                fprintf(' 通过\n');
            else
                fprintf(' 失败\n');
            end
            
            % 测试3: 约束强制
            fprintf('测试3: 约束强制...');
            X_corrected = ManifoldOperations.enforce_constraints(X_invalid);
            valid3 = ManifoldOperations.is_valid_state(X_corrected);
            
            if valid3
                fprintf(' 通过\n');
            else
                fprintf(' 失败\n');
            end
            
            % 测试4: 雅可比矩阵维度
            fprintf('测试4: 雅可比矩阵维度...');
            xi_test = randn(ManifoldOperations.TANGENT_DIM, 1) * 0.01;
            J = ManifoldOperations.jacobian_boxplus(X1, xi_test);
            
            correct_size = (size(J, 1) == ManifoldOperations.STATE_DIM) && ...
                          (size(J, 2) == ManifoldOperations.TANGENT_DIM);
            
            if correct_size
                fprintf(' 通过 (大小: %dx%d)\n', size(J, 1), size(J, 2));
            else
                fprintf(' 失败 (大小: %dx%d)\n', size(J, 1), size(J, 2));
            end
            
            fprintf('=== 测试完成 ===\n\n');
        end
        
    end
    
end

function str = logical_to_string(val)
% 将逻辑值转换为字符串
if val
    str = '有效';
else
    str = '无效';
end
end