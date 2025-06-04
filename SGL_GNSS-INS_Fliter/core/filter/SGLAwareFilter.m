classdef SGLAwareFilter < handle
% SGLAWAREFILTER SGL-aware自适应滤波器核心类
%
% 
% 包括四种运行模式：ACTIVE, OBSERVE, FROZEN, REACTIVATE
%
% 核心功能:
% - 复合流形上的几何随机微分方程建模
% - Fisher信息矩阵的几何增强计算
% - SGL检测与四模式智能切换
% - 自适应噪声参数估计
% - 几何约束保持的状态估计
%
% 作者: Y
% 日期: 2025.6

    properties (Access = private)
        % 配置参数
        config
        
        % 当前状态
        X_current           % 17x1 当前状态估计
        P_current           % 15x15 当前协方差矩阵(切空间)
        
        % 噪声参数
        theta_current       % 当前噪声参数向量
        Q_current           % 过程噪声协方差矩阵
        R_current           % 观测噪声协方差矩阵
        
        % SGL相关
        fisher_info         % 当前Fisher信息矩阵
        sgl_detector        % SGL检测器
        mode_switcher       % 模式切换器
        current_mode        % 当前运行模式 (1=ACTIVE, 2=OBSERVE, 3=FROZEN, 4=REACTIVATE)
        
        % 历史数据 (用于自适应)
        innovation_history  % 新息序列历史
        residual_history    % 残差历史
        info_history        % 信息历史
        
        % 性能统计
        compute_times       % 计算时间统计
        mode_counts         % 模式计数
        
        % 几何模型
        geometric_predictor % 几何预测器
        geometric_updater   % 几何更新器
        
        % 时间管理
        current_time        % 当前时间
        last_update_time    % 上次更新时间
        
        % 调试信息
        debug_info          % 调试信息结构
    end
    
    properties (Constant)
        % 运行模式常数
        MODE_ACTIVE = 1;
        MODE_OBSERVE = 2;
        MODE_FROZEN = 3;
        MODE_REACTIVATE = 4;
        
        % 模式名称
        MODE_NAMES = {'ACTIVE', 'OBSERVE', 'FROZEN', 'REACTIVATE'};
    end
    
    methods (Access = public)
        
        function obj = SGLAwareFilter(config)
            % 构造函数
            %
            % 输入:
            %   config - 配置结构体
            
            obj.config = config;
            
            % 初始化各模块
            obj.init_sgl_detector();
            obj.init_mode_switcher();
            obj.init_geometric_modules();
            obj.init_noise_parameters();
            obj.init_history_buffers();
            obj.init_performance_counters();
            
            % 设置初始模式
            obj.current_mode = obj.MODE_ACTIVE;
            obj.current_time = 0;
            obj.last_update_time = 0;
            
            if obj.config.debug.verbose
                fprintf('SGL-aware滤波器初始化完成\n');
            end
        end
        
        function initialize(obj, X0, P0)
            % 初始化滤波器状态
            %
            % 输入:
            %   X0 - 17x1 初始状态向量
            %   P0 - 15x15 初始协方差矩阵(切空间) 或 17x17(状态空间)
            
            % 验证输入
            if ~ManifoldOperations.is_valid_state(X0)
                warning('初始状态无效，正在修正...');
                X0 = ManifoldOperations.enforce_constraints(X0);
            end
            
            obj.X_current = X0;
            
            % 处理协方差矩阵维度
            if size(P0, 1) == 17 && size(P0, 2) == 17
                % 状态空间协方差，需要投影到切空间
                obj.P_current = obj.project_covariance_to_tangent(P0);
            else
                obj.P_current = P0;
            end
            
            % 初始化Fisher信息矩阵
            obj.fisher_info = obj.compute_initial_fisher_info();
            
            if obj.config.debug.verbose
                fprintf('滤波器状态初始化完成\n');
                ManifoldOperations.display_state(obj.X_current, '初始状态');
            end
        end
        
        function predict(obj, imu_data)
            % 预测步骤 - 基于IMU数据进行状态预测
            %
            % 输入:
            %   imu_data - 结构体或数组，包含IMU测量数据
            %     .timestamp - 时间戳
            %     .acc - 3x1 加速度计测量 [m/s²]
            %     .gyro - 3x1 陀螺仪测量 [rad/s]
            
            tic; % 开始计时
            
            % 解析IMU数据
            if isstruct(imu_data)
                timestamp = imu_data.timestamp;
                acc_meas = imu_data.acc(:);
                gyro_meas = imu_data.gyro(:);
            else
                % 假设格式: [time, ax, ay, az, gx, gy, gz]
                timestamp = imu_data(1);
                acc_meas = imu_data(2:4)';
                gyro_meas = imu_data(5:7)';
            end
            
            % 计算时间间隔
            if obj.current_time > 0
                dt = timestamp - obj.current_time;
            else
                dt = obj.config.dt.imu;
            end
            
            % 几何预测 (连续时间几何随机微分方程)
            [X_pred, P_pred] = obj.geometric_predictor.predict(...
                obj.X_current, obj.P_current, acc_meas, gyro_meas, dt, obj.Q_current);
            
            % 确保几何约束
            X_pred = ManifoldOperations.enforce_constraints(X_pred);
            
            % 更新状态
            obj.X_current = X_pred;
            obj.P_current = P_pred;
            obj.current_time = timestamp;
            
            % 记录性能
            obj.compute_times.predict = [obj.compute_times.predict, toc];
            
            if obj.config.debug.verbose && mod(length(obj.compute_times.predict), 1000) == 0
                fprintf('预测步骤 %d 完成，平均耗时: %.3f ms\n', ...
                        length(obj.compute_times.predict), ...
                        mean(obj.compute_times.predict) * 1000);
            end
        end
        
        function update(obj, gnss_data)
            % 更新步骤 - 基于GNSS观测进行状态更新
            %
            % 输入:
            %   gnss_data - 结构体或数组，包含GNSS观测数据
            %     .timestamp - 时间戳
            %     .pos - 3x1 位置观测 [m]
            %     .vel - 3x1 速度观测 [m/s] (可选)
            %     .std_pos - 3x1 位置标准差 [m] (可选)
            %     .std_vel - 3x1 速度标准差 [m/s] (可选)
            
            tic; % 开始计时
            
            % 解析GNSS数据
            obs_data = obj.parse_gnss_data(gnss_data);
            
            % 更新时间
            obj.last_update_time = obj.current_time;
            
            % 几何更新
            [X_updated, P_updated, innovation, S] = obj.geometric_updater.update(...
                obj.X_current, obj.P_current, obs_data, obj.R_current);
            
            % 确保几何约束
            X_updated = ManifoldOperations.enforce_constraints(X_updated);
            
            % 更新Fisher信息矩阵
            obj.update_fisher_information(innovation, S);
            
            % SGL检测和模式控制
            obj.detect_sgl_and_switch_mode(innovation, S);
            
            % 自适应噪声参数更新 (根据当前模式)
            if obj.current_mode ~= obj.MODE_FROZEN
                obj.update_noise_parameters(innovation, S);
            end
            
            % 更新状态
            obj.X_current = X_updated;
            obj.P_current = P_updated;
            
            % 记录历史数据
            obj.record_history(innovation, S);
            
            % 记录性能
            obj.compute_times.update = [obj.compute_times.update, toc];
            obj.mode_counts(obj.current_mode) = obj.mode_counts(obj.current_mode) + 1;
            
            if obj.config.debug.verbose
                fprintf('时间 %.2f: 更新完成, 模式=%s, 新息范数=%.4f\n', ...
                        obj.current_time, obj.MODE_NAMES{obj.current_mode}, norm(innovation));
            end
        end
        
        function state = get_current_state(obj)
            % 获取当前状态
            %
            % 输出:
            %   state - 结构体，包含当前状态信息
            
            [pos, vel, att_euler, ba, bg, clk] = ...
                ManifoldOperations.extract_state_components(obj.X_current);
            
            state = struct();
            state.timestamp = obj.current_time;
            state.X = obj.X_current;
            state.position = pos;
            state.velocity = vel;
            state.attitude_euler = att_euler;
            state.bias_acc = ba;
            state.bias_gyro = bg;
            state.clock = clk;
            state.mode = obj.current_mode;
            state.mode_name = obj.MODE_NAMES{obj.current_mode};
        end
        
        function cov = get_current_covariance(obj)
            % 获取当前协方差矩阵
            %
            % 输出:
            %   cov - 15x15 协方差矩阵(切空间)
            
            cov = obj.P_current;
        end
        
        function fisher_info = get_fisher_info(obj)
            % 获取当前Fisher信息矩阵
            %
            % 输出:
            %   fisher_info - Fisher信息矩阵
            
            fisher_info = obj.fisher_info;
        end
        
        function noise_params = get_noise_params(obj)
            % 获取当前噪声参数
            %
            % 输出:
            %   noise_params - 噪声参数结构体
            
            noise_params = struct();
            noise_params.theta = obj.theta_current;
            noise_params.Q = obj.Q_current;
            noise_params.R = obj.R_current;
        end
        
        function mode = get_current_mode(obj)
            % 获取当前运行模式
            %
            % 输出:
            %   mode - 当前模式编号
            
            mode = obj.current_mode;
        end
        
        function metrics = get_performance_metrics(obj)
            % 获取性能指标
            %
            % 输出:
            %   metrics - 10x1 性能指标向量
            
            metrics = zeros(10, 1);
            
            % 计算负荷 (基于当前模式)
            load_factors = [1.0, 0.3, 0.05, 1.2]; % ACTIVE, OBSERVE, FROZEN, REACTIVATE
            metrics(1) = load_factors(obj.current_mode) * 100; % 百分比
            
            % 平均计算时间
            if ~isempty(obj.compute_times.predict)
                metrics(2) = mean(obj.compute_times.predict) * 1000; % ms
            end
            if ~isempty(obj.compute_times.update)
                metrics(3) = mean(obj.compute_times.update) * 1000; % ms
            end
            
            % Fisher信息统计
            if ~isempty(obj.fisher_info)
                metrics(4) = trace(obj.fisher_info);
                metrics(5) = min(eig(obj.fisher_info));
                metrics(6) = cond(obj.fisher_info);
            end
            
            % 模式分布
            total_updates = sum(obj.mode_counts);
            if total_updates > 0
                metrics(7) = obj.mode_counts(obj.MODE_FROZEN) / total_updates * 100; % FROZEN百分比
            end
            
            % 新息统计
            if ~isempty(obj.innovation_history)
                recent_innovations = obj.innovation_history(max(1, end-99):end, :);
                metrics(8) = mean(sqrt(sum(recent_innovations.^2, 2))); % 平均新息范数
            end
            
            % 协方差迹
            metrics(9) = trace(obj.P_current);
            
            % SGL饱和度指标
            if ~isempty(obj.fisher_info)
                lambda_min = min(eig(obj.fisher_info));
                theoretical_max = obj.config.computed.sgl_saturation_bound;
                metrics(10) = lambda_min / theoretical_max * 100; % 饱和度百分比
            end
        end
        
    end
    
    methods (Access = private)
        
        function init_sgl_detector(obj)
            % 初始化SGL检测器
            
            obj.sgl_detector = SGLDetector(obj.config.sgl);
        end
        
        function init_mode_switcher(obj)
            % 初始化模式切换器
            
            obj.mode_switcher = ModeSwitcher(obj.config.modes);
        end
        
        function init_geometric_modules(obj)
            % 初始化几何模块
            
            obj.geometric_predictor = GeometricPredictor(obj.config);
            obj.geometric_updater = GeometricUpdater(obj.config);
        end
        
        function init_noise_parameters(obj)
            % 初始化噪声参数
            
            % 过程噪声协方差矩阵
            obj.Q_current = obj.config.computed.Q_discrete;
            
            % 观测噪声协方差矩阵
            obj.R_current = obj.config.computed.R_gnss;
            
            % 噪声参数向量 (用于自适应)
            obj.theta_current = obj.pack_noise_parameters();
        end
        
        function init_history_buffers(obj)
            % 初始化历史数据缓冲区
            
            max_history = 1000; % 最大历史长度
            
            obj.innovation_history = [];
            obj.residual_history = [];
            obj.info_history = [];
        end
        
        function init_performance_counters(obj)
            % 初始化性能计数器
            
            obj.compute_times = struct();
            obj.compute_times.predict = [];
            obj.compute_times.update = [];
            
            obj.mode_counts = zeros(4, 1); % 四种模式的计数
        end
        
        function theta = pack_noise_parameters(obj)
            % 打包噪声参数为向量
            
            % 提取协方差矩阵的独立元素
            Q_vec = obj.extract_covariance_params(obj.Q_current);
            R_vec = obj.extract_covariance_params(obj.R_current);
            
            theta = [Q_vec; R_vec];
        end
        
        function params = extract_covariance_params(obj, C)
            % 提取协方差矩阵的独立参数
            %
            % 对于对角矩阵，提取对角元素
            % 对于一般矩阵，提取上三角元素
            
            if obj.is_diagonal_matrix(C)
                params = diag(C);
            else
                % 提取上三角元素 (包括对角线)
                n = size(C, 1);
                params = C(triu(true(n)));
            end
        end
        
        function is_diag = is_diagonal_matrix(obj, M)
            % 检查矩阵是否为对角矩阵
            
            is_diag = norm(M - diag(diag(M)), 'fro') < 1e-10;
        end
        
        function obs_data = parse_gnss_data(obj, gnss_data)
            % 解析GNSS观测数据
            
            if isstruct(gnss_data)
                obs_data = gnss_data;
            else
                % 假设格式: [time, x, y, z, vx, vy, vz, ...]
                obs_data = struct();
                obs_data.timestamp = gnss_data(1);
                obs_data.pos = gnss_data(2:4)';
                if length(gnss_data) >= 7
                    obs_data.vel = gnss_data(5:7)';
                else
                    obs_data.vel = [];
                end
                if length(gnss_data) >= 10
                    obs_data.std_pos = gnss_data(8:10)';
                else
                    obs_data.std_pos = [];
                end
            end
        end
        
        function update_fisher_information(obj, innovation, S)
            % 更新Fisher信息矩阵
            %
            % 实现文章中的几何增强Fisher信息矩阵计算
            
            % 观测信息贡献
            H = obj.geometric_updater.get_last_observation_matrix();
            I_obs = H' * (S \ H); % H^T * S^(-1) * H
            
            % 流形信息贡献 (基于当前参数的几何结构)
            I_manifold = obj.compute_manifold_information();
            
            % 几何增强Fisher信息矩阵 (文章公式)
            I_geo = I_obs + I_manifold;
            
            % 累积更新 (考虑遗忘因子)
            forgetting_factor = obj.config.adaptive.forgetting_factor;
            if isempty(obj.fisher_info)
                obj.fisher_info = I_geo;
            else
                obj.fisher_info = forgetting_factor * obj.fisher_info + I_geo;
            end
        end
        
        function I_manifold = compute_manifold_information(obj)
            % 计算流形信息贡献
            %
            % 基于参数空间的几何结构计算信息贡献
            
            % 这里简化实现，实际应根据文章中的流形几何计算
            param_dim = length(obj.theta_current);
            I_manifold = 0.01 * eye(param_dim); % 基础流形信息
            
            % 根据当前噪声参数调整
            % 对于正定矩阵参数，使用Fisher-Rao度量
            % 这里需要更详细的实现，暂时简化
        end
        
        function detect_sgl_and_switch_mode(obj, innovation, S)
            % SGL检测和模式切换
            
            % SGL检测
            sgl_detected = obj.sgl_detector.detect_sgl(obj.fisher_info, innovation, S);
            
            % 环境变化检测
            chi2_stat = innovation' * (S \ innovation);
            env_changed = chi2_stat > obj.config.modes.switch.chi2_threshold;
            
            % 模式决策
            new_mode = obj.mode_switcher.decide_mode(...
                obj.current_mode, sgl_detected, env_changed, norm(innovation));
            
            % 模式切换
            if new_mode ~= obj.current_mode
                if obj.config.debug.verbose
                    fprintf('模式切换: %s -> %s\n', ...
                            obj.MODE_NAMES{obj.current_mode}, obj.MODE_NAMES{new_mode});
                end
                obj.current_mode = new_mode;
                
                % 模式切换时的特殊处理
                obj.handle_mode_transition(new_mode);
            end
        end
        
        function handle_mode_transition(obj, new_mode)
            % 处理模式切换
            
            switch new_mode
                case obj.MODE_ACTIVE
                    % 激活模式：重置某些参数
                    obj.sgl_detector.reset_counters();
                    
                case obj.MODE_OBSERVE
                    % 观察模式：降低更新频率
                    
                case obj.MODE_FROZEN
                    % 冻结模式：停止参数更新
                    if obj.config.debug.verbose
                        fprintf('进入FROZEN模式，参数冻结\n');
                    end
                    
                case obj.MODE_REACTIVATE
                    % 重激活模式：重新启动自适应
                    obj.sgl_detector.reset_counters();
                    if obj.config.debug.verbose
                        fprintf('环境变化检测，重新激活自适应\n');
                    end
            end
        end
        
        function update_noise_parameters(obj, innovation, S)
            % 自适应噪声参数更新
            %
            % 基于自然梯度方法更新噪声参数
            
            if obj.current_mode == obj.MODE_FROZEN
                return; % 冻结模式不更新参数
            end
            
            % 计算步长 (基于当前模式和Fisher信息)
            alpha = obj.compute_adaptive_step_size();
            
            % 计算梯度
            grad = obj.compute_parameter_gradient(innovation, S);
            
            % 自然梯度更新 (简化实现)
            if ~isempty(obj.fisher_info) && cond(obj.fisher_info) < 1e12
                natural_grad = obj.fisher_info \ grad;
            else
                natural_grad = grad; % 回退到普通梯度
            end
            
            % 参数更新
            obj.theta_current = obj.theta_current - alpha * natural_grad;
            
            % 解包参数并更新协方差矩阵
            obj.unpack_noise_parameters();
            
            % 确保协方差矩阵正定性
            obj.enforce_covariance_constraints();
        end
        
        function alpha = compute_adaptive_step_size(obj)
            % 计算自适应步长
            %
            % 实现文章中的步长计算公式
            
            % 基础步长
            mode_factors = [1.0, 0.3, 0, 1.2]; % 对应四种模式
            alpha_base = obj.config.modes.active.alpha_base;
            
            % 模式调制
            alpha = alpha_base * mode_factors(obj.current_mode);
            
            % Fisher信息调制
            if ~isempty(obj.fisher_info)
                fisher_trace = trace(obj.fisher_info);
                C_step = obj.config.adaptive.natural_gradient.C_step;
                alpha = alpha * min(1, C_step / sqrt(fisher_trace));
            end
        end
        
        function grad = compute_parameter_gradient(obj, innovation, S)
            % 计算参数梯度
            %
            % 这里简化实现，实际应根据具体的似然函数计算
            
            param_dim = length(obj.theta_current);
            grad = zeros(param_dim, 1);
            
            % 基于新息的梯度估计 (简化)
            chi2 = innovation' * (S \ innovation);
            obs_dim = length(innovation);
            
            % 观测噪声参数梯度
            R_dim = size(obj.R_current, 1);
            R_grad_scale = (chi2 - obs_dim) / obs_dim;
            
            % 这里需要更精确的梯度计算，暂时简化
            grad(end-R_dim+1:end) = R_grad_scale * ones(R_dim, 1) * 0.01;
        end
        
        function unpack_noise_parameters(obj)
            % 解包噪声参数
            
            % 这里需要根据参数向量重构协方差矩阵
            % 简化实现：假设参数直接对应协方差矩阵的对角元素
            
            Q_dim = size(obj.Q_current, 1);
            R_dim = size(obj.R_current, 1);
            
            if length(obj.theta_current) >= Q_dim + R_dim
                Q_params = obj.theta_current(1:Q_dim);
                R_params = obj.theta_current(Q_dim+1:Q_dim+R_dim);
                
                % 确保正值
                Q_params = max(Q_params, 1e-6);
                R_params = max(R_params, 1e-6);
                
                obj.Q_current = diag(Q_params);
                obj.R_current = diag(R_params);
            end
        end
        
        function enforce_covariance_constraints(obj)
            % 强制协方差矩阵约束
            
            % 确保正定性
            [V, D] = eig(obj.Q_current);
            D = diag(max(diag(D), 1e-8));
            obj.Q_current = V * D * V';
            
            [V, D] = eig(obj.R_current);
            D = diag(max(diag(D), 1e-8));
            obj.R_current = V * D * V';
            
            % 确保数值稳定性
            obj.Q_current = (obj.Q_current + obj.Q_current') / 2;
            obj.R_current = (obj.R_current + obj.R_current') / 2;
        end
        
        function record_history(obj, innovation, S)
            % 记录历史数据
            
            max_history = 1000;
            
            % 新息历史
            obj.innovation_history = [obj.innovation_history; innovation'];
            if size(obj.innovation_history, 1) > max_history
                obj.innovation_history = obj.innovation_history(end-max_history+1:end, :);
            end
            
            % 残差统计历史
            chi2_stat = innovation' * (S \ innovation);
            obj.residual_history = [obj.residual_history; chi2_stat];
            if length(obj.residual_history) > max_history
                obj.residual_history = obj.residual_history(end-max_history+1:end);
            end
            
            % Fisher信息历史
            if ~isempty(obj.fisher_info)
                fisher_trace = trace(obj.fisher_info);
                obj.info_history = [obj.info_history; fisher_trace];
                if length(obj.info_history) > max_history
                    obj.info_history = obj.info_history(end-max_history+1:end);
                end
            end
        end
        
        function P_tangent = project_covariance_to_tangent(obj, P_state)
            % 将状态空间协方差投影到切空间
            %
            % 输入:
            %   P_state - 17x17 状态空间协方差矩阵
            %
            % 输出:
            %   P_tangent - 15x15 切空间协方差矩阵
            
            % 提取切空间相关的协方差块
            % 位置、速度、姿态(切空间)、IMU偏差
            idx_tangent = [1:6, 7:9, 10:15]; % 排除时钟参数
            P_tangent = P_state(idx_tangent, idx_tangent);
        end
        
        function I_init = compute_initial_fisher_info(obj)
            % 计算初始Fisher信息矩阵
            
            param_dim = length(obj.theta_current);
            
            % 基于初始协方差的Fisher信息
            I_init = obj.config.adaptive.min_samples * eye(param_dim);
            
            % 添加先验信息
            prior_strength = 0.1;
            I_init = I_init + prior_strength * eye(param_dim);
        end
        
    end
    
    methods (Static)
        
        function run_tests()
            % 运行SGL-aware滤波器的单元测试
            
            fprintf('=== SGL-aware滤波器单元测试 ===\n');
            
            % 创建测试配置
            config = config_system();
            
            % 测试1: 滤波器初始化
            fprintf('测试1: 滤波器初始化...');
            try
                filter = SGLAwareFilter(config);
                fprintf(' 通过\n');
            catch ME
                fprintf(' 失败: %s\n', ME.message);
                return;
            end
            
            % 测试2: 状态初始化
            fprintf('测试2: 状态初始化...');
            try
                X0 = ManifoldOperations.create_initial_state(...
                    [0;0;0], [0;0;0], [0;0;0], [0;0;0], [0;0;0], [0;0]);
                P0 = config.computed.P0;
                filter.initialize(X0, P0);
                fprintf(' 通过\n');
            catch ME
                fprintf(' 失败: %s\n', ME.message);
                return;
            end
            
            % 测试3: 预测步骤
            fprintf('测试3: 预测步骤...');
            try
                imu_data = [0.1, 0, 0, -9.8, 0, 0, 0]; % [time, acc, gyro]
                filter.predict(imu_data);
                fprintf(' 通过\n');
            catch ME
                fprintf(' 失败: %s\n', ME.message);
                return;
            end
            
            % 测试4: 更新步骤
            fprintf('测试4: 更新步骤...');
            try
                gnss_data = [0.1, 1, 1, 0, 0.1, 0.1, 0]; % [time, pos, vel]
                filter.update(gnss_data);
                fprintf(' 通过\n');
            catch ME
                fprintf(' 失败: %s\n', ME.message);
                return;
            end
            
            % 测试5: 状态获取
            fprintf('测试5: 状态获取...');
            try
                state = filter.get_current_state();
                cov = filter.get_current_covariance();
                metrics = filter.get_performance_metrics();
                fprintf(' 通过\n');
            catch ME
                fprintf(' 失败: %s\n', ME.message);
                return;
            end
            
            fprintf('=== 所有测试通过 ===\n\n');
        end
        
    end
    
end