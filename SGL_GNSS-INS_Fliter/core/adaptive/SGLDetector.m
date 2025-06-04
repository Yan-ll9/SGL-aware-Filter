classdef SGLDetector < handle
% SGLDETECTOR 标准几何极限(SGL)检测器
%
% SGL理论检测机制，包括四个检测准则：
% 1. 谱饱和准则 - Fisher信息矩阵特征值分布均匀性
% 2. 条件数稳定准则 - 参数相关性收敛  
% 3. 信息增量衰减准则 - 信息增长率衰减
% 4. 参数收敛准则 - 参数估计稳定性
%
%
% 作者: Y
% 日期: 2025.6

    properties (Access = private)
        % 配置参数
        config
        
        % 检测阈值
        tau1            % 谱饱和准则阈值
        tau2            % 条件数稳定准则阈值
        tau3            % 信息增量衰减准则阈值
        tau4            % 参数收敛准则阈值
        
        % 历史数据窗口
        window_length   % 观察窗口长度
        confirm_count   % 确认阈值
        
        % 内部状态
        fisher_history          % Fisher信息矩阵历史
        innovation_history      % 新息历史
        parameter_history       % 参数历史
        
        % 检测状态
        sgl_count              % 连续检测到SGL的次数
        last_detection_time    % 上次检测时间
        detection_flags        % 各准则检测标志
        
        % 统计信息
        detection_stats        % 检测统计
        performance_metrics    % 性能指标
        
        % 理论参数
        C_sat              % 饱和常数
        alpha_theory       % 几何约束参数α
        epsilon_theory     % 几何约束参数ε
        
        % 调试信息
        debug_info         % 调试信息记录
    end
    
    methods (Access = public)
        
        function obj = SGLDetector(config)
            % 构造函数
            %
            % 输入:
            %   config - SGL配置结构体
            
            obj.config = config;
            
            % 初始化检测阈值
            obj.tau1 = config.detection.tau1;
            obj.tau2 = config.detection.tau2;
            obj.tau3 = config.detection.tau3;
            obj.tau4 = config.detection.tau4;
            
            % 初始化窗口参数
            obj.window_length = config.detection.window_length;
            obj.confirm_count = config.detection.confirm_count;
            
            % 初始化理论参数
            obj.C_sat = config.theory.C_sat;
            obj.alpha_theory = config.theory.alpha;
            obj.epsilon_theory = config.theory.epsilon;
            
            % 初始化内部状态
            obj.reset_state();
            
            % 初始化统计信息
            obj.init_statistics();
        end
        
        function sgl_detected = detect_sgl(obj, fisher_info, innovation, S)
            % 主检测函数 - 基于四准则综合判定
            %
            % 输入:
            %   fisher_info - 当前Fisher信息矩阵
            %   innovation - 新息向量
            %   S - 新息协方差矩阵
            %
            % 输出:
            %   sgl_detected - 逻辑值，是否检测到SGL
            
            % 更新历史数据
            obj.update_history(fisher_info, innovation);
            
            % 计算四个检测准则
            C1 = obj.compute_spectral_saturation_criterion(fisher_info);
            C2 = obj.compute_condition_number_criterion(fisher_info);
            C3 = obj.compute_information_increment_criterion();
            C4 = obj.compute_parameter_convergence_criterion();
            
            % 记录检测标志
            obj.detection_flags = [C1, C2, C3, C4];
            
            % 综合SGL检测逻辑
            % 核心条件: C3 AND C4 (信息增量衰减 且 参数收敛)
            % 辅助条件: C1 OR C2 (谱饱和 或 条件数稳定)
            core_condition = C3 && C4;
            auxiliary_condition = C1 || C2;
            
            current_detection = core_condition && auxiliary_condition;
            
            % 连续检测确认机制
            if current_detection
                obj.sgl_count = obj.sgl_count + 1;
            else
                obj.sgl_count = 0;
            end
            
            % 最终SGL检测结果
            sgl_detected = obj.sgl_count >= obj.confirm_count;
            
            % 更新统计信息
            obj.update_detection_statistics(sgl_detected, obj.detection_flags);
            
            % 调试信息
            if obj.config.debug.verbose && sgl_detected
                obj.log_detection_event(fisher_info, innovation);
            end
        end
        
        function C1 = compute_spectral_saturation_criterion(obj, fisher_info)
            % 准则1: 谱饱和准则
            %
            % 数学表达式: C1 = λ_min(I_geo) / (tr(I_geo)/dim(θ)) ≥ τ1
            % 理论依据: 特征值谱均匀分布表明信息积累达到平衡
            %
            % 输入:
            %   fisher_info - Fisher信息矩阵
            %
            % 输出:
            %   C1 - 逻辑值，是否满足谱饱和准则
            
            if isempty(fisher_info) || any(~isfinite(fisher_info(:)))
                C1 = false;
                return;
            end
            
            % 计算特征值
            eigenvals = eig(fisher_info);
            
            % 确保特征值为实数且非负
            eigenvals = real(eigenvals);
            eigenvals = max(eigenvals, 0);
            
            if min(eigenvals) < 1e-12
                C1 = false;
                return;
            end
            
            % 计算谱饱和指标
            lambda_min = min(eigenvals);
            trace_fisher = trace(fisher_info);
            dim_theta = size(fisher_info, 1);
            
            if trace_fisher < 1e-12 || dim_theta == 0
                C1 = false;
                return;
            end
            
            spectral_ratio = lambda_min / (trace_fisher / dim_theta);
            
            C1 = spectral_ratio >= obj.tau1;
            
            % 记录调试信息
            obj.debug_info.last_spectral_ratio = spectral_ratio;
        end
        
        function C2 = compute_condition_number_criterion(obj, fisher_info)
            % 准则2: 条件数稳定准则
            %
            % 数学表达式: C2 = cond(I_geo) ≤ τ2
            % 理论依据: 参数相关性收敛，条件数趋于稳定
            %
            % 输入:
            %   fisher_info - Fisher信息矩阵
            %
            % 输出:
            %   C2 - 逻辑值，是否满足条件数稳定准则
            
            if isempty(fisher_info) || any(~isfinite(fisher_info(:)))
                C2 = false;
                return;
            end
            
            % 计算条件数
            condition_number = cond(fisher_info);
            
            % 处理病态矩阵
            if ~isfinite(condition_number) || condition_number > 1e15
                C2 = false;
                return;
            end
            
            C2 = condition_number <= obj.tau2;
            
            % 记录调试信息
            obj.debug_info.last_condition_number = condition_number;
        end
        
        function C3 = compute_information_increment_criterion(obj)
            % 准则3: 信息增量衰减准则
            %
            % 数学表达式: C3 = ||I_k - I_{k-1}||_F / ||I_{k-1}||_F ≤ τ3
            % 理论依据: 信息增长率衰减表明接近饱和状态
            %
            % 输出:
            %   C3 - 逻辑值，是否满足信息增量衰减准则
            
            if length(obj.fisher_history) < 2
                C3 = false;
                return;
            end
            
            % 获取最近的两个Fisher信息矩阵
            I_current = obj.fisher_history{end};
            I_previous = obj.fisher_history{end-1};
            
            if isempty(I_current) || isempty(I_previous) || ...
               any(~isfinite(I_current(:))) || any(~isfinite(I_previous(:)))
                C3 = false;
                return;
            end
            
            % 计算Frobenius范数的相对变化
            norm_previous = norm(I_previous, 'fro');
            
            if norm_previous < 1e-12
                C3 = false;
                return;
            end
            
            norm_increment = norm(I_current - I_previous, 'fro');
            relative_increment = norm_increment / norm_previous;
            
            C3 = relative_increment <= obj.tau3;
            
            % 记录调试信息
            obj.debug_info.last_info_increment = relative_increment;
        end
        
        function C4 = compute_parameter_convergence_criterion(obj)
            % 准则4: 参数收敛准则
            %
            % 数学表达式: C4 = ||θ_k - θ_{k-L}||_Θ / ||θ_{k-L}||_Θ ≤ τ4
            % 理论依据: 参数估计稳定表明收敛到真值邻域
            %
            % 输出:
            %   C4 - 逻辑值，是否满足参数收敛准则
            
            if length(obj.parameter_history) < obj.window_length
                C4 = false;
                return;
            end
            
            % 获取当前参数和L步前的参数
            theta_current = obj.parameter_history{end};
            theta_past = obj.parameter_history{end - obj.window_length + 1};
            
            if isempty(theta_current) || isempty(theta_past) || ...
               any(~isfinite(theta_current)) || any(~isfinite(theta_past))
                C4 = false;
                return;
            end
            
            % 计算参数空间Θ上的距离（这里简化为欧几里得距离）
            norm_past = norm(theta_past);
            
            if norm_past < 1e-12
                C4 = false;
                return;
            end
            
            param_change = norm(theta_current - theta_past) / norm_past;
            
            C4 = param_change <= obj.tau4;
            
            % 记录调试信息
            obj.debug_info.last_param_change = param_change;
        end
        
        function theoretical_bound = compute_theoretical_sgl_bound(obj, n)
            % 计算理论SGL界限
            %
            % 基于文章中的Fisher信息饱和定理计算理论上界
            %
            % 输入:
            %   n - 观测次数
            %
            % 输出:
            %   theoretical_bound - 理论SGL界限
            
            if nargin < 2
                n = length(obj.fisher_history);
            end
            
            % 实现文章公式: λ_min(I_geo^(n)) ≤ C_sat · n^α / (1 + ε·n^α)
            theoretical_bound = obj.C_sat * (n^obj.alpha_theory) / ...
                               (1 + obj.epsilon_theory * (n^obj.alpha_theory));
        end
        
        function saturation_ratio = compute_saturation_ratio(obj, fisher_info)
            % 计算当前饱和度比例
            %
            % 输入:
            %   fisher_info - 当前Fisher信息矩阵
            %
            % 输出:
            %   saturation_ratio - 饱和度比例 [0, 1]
            
            if isempty(fisher_info)
                saturation_ratio = 0;
                return;
            end
            
            % 当前信息量
            current_info = min(eig(fisher_info));
            
            % 理论最大值
            n = length(obj.fisher_history);
            max_theoretical = obj.compute_theoretical_sgl_bound(n);
            
            % 饱和度比例
            if max_theoretical > 1e-12
                saturation_ratio = min(current_info / max_theoretical, 1.0);
            else
                saturation_ratio = 0;
            end
        end
        
        function reset_counters(obj)
            % 重置检测计数器
            
            obj.sgl_count = 0;
            obj.last_detection_time = 0;
        end
        
        function reset_state(obj)
            % 重置检测器状态
            
            obj.fisher_history = {};
            obj.innovation_history = [];
            obj.parameter_history = {};
            
            obj.sgl_count = 0;
            obj.last_detection_time = 0;
            obj.detection_flags = [false, false, false, false];
            
            obj.debug_info = struct();
            obj.debug_info.last_spectral_ratio = 0;
            obj.debug_info.last_condition_number = inf;
            obj.debug_info.last_info_increment = inf;
            obj.debug_info.last_param_change = inf;
        end
        
        function stats = get_detection_statistics(obj)
            % 获取检测统计信息
            %
            % 输出:
            %   stats - 检测统计结构体
            
            stats = obj.detection_stats;
            
            % 添加当前状态
            stats.current_sgl_count = obj.sgl_count;
            stats.detection_flags = obj.detection_flags;
            stats.debug_info = obj.debug_info;
        end
        
        function display_status(obj)
            % 显示检测器状态（调试用）
            
            fprintf('\n=== SGL检测器状态 ===\n');
            
            fprintf('检测阈值: τ1=%.3f, τ2=%.1f, τ3=%.3f, τ4=%.3f\n', ...
                    obj.tau1, obj.tau2, obj.tau3, obj.tau4);
            
            fprintf('理论参数: C_sat=%.2f, α=%.2f, ε=%.2f\n', ...
                    obj.C_sat, obj.alpha_theory, obj.epsilon_theory);
            
            fprintf('当前状态:\n');
            fprintf('  连续SGL检测次数: %d/%d\n', obj.sgl_count, obj.confirm_count);
            
            if ~isempty(obj.detection_flags)
                flag_names = {'谱饱和', '条件数稳定', '信息增量衰减', '参数收敛'};
                for i = 1:4
                    status = obj.detection_flags(i);
                    fprintf('  %s: %s\n', flag_names{i}, logical_to_string(status));
                end
            end
            
            fprintf('调试信息:\n');
            if isfield(obj.debug_info, 'last_spectral_ratio')
                fprintf('  谱比例: %.4f (阈值: %.3f)\n', ...
                        obj.debug_info.last_spectral_ratio, obj.tau1);
            end
            if isfield(obj.debug_info, 'last_condition_number')
                fprintf('  条件数: %.2e (阈值: %.1f)\n', ...
                        obj.debug_info.last_condition_number, obj.tau2);
            end
            if isfield(obj.debug_info, 'last_info_increment')
                fprintf('  信息增量: %.4f (阈值: %.3f)\n', ...
                        obj.debug_info.last_info_increment, obj.tau3);
            end
            if isfield(obj.debug_info, 'last_param_change')
                fprintf('  参数变化: %.4f (阈值: %.3f)\n', ...
                        obj.debug_info.last_param_change, obj.tau4);
            end
            
            fprintf('====================\n\n');
        end
        
    end
    
    methods (Access = private)
        
        function update_history(obj, fisher_info, innovation)
            % 更新历史数据
            
            % 更新Fisher信息历史
            obj.fisher_history{end+1} = fisher_info;
            if length(obj.fisher_history) > obj.window_length
                obj.fisher_history = obj.fisher_history(end-obj.window_length+1:end);
            end
            
            % 更新新息历史
            obj.innovation_history = [obj.innovation_history; innovation'];
            if size(obj.innovation_history, 1) > obj.window_length
                obj.innovation_history = obj.innovation_history(end-obj.window_length+1:end, :);
            end
        end
        
        function init_statistics(obj)
            % 初始化统计信息
            
            obj.detection_stats = struct();
            obj.detection_stats.total_detections = 0;
            obj.detection_stats.false_positive_rate = 0;
            obj.detection_stats.detection_delay = [];
            obj.detection_stats.criterion_success_rate = zeros(4, 1);
            
            obj.performance_metrics = struct();
            obj.performance_metrics.average_detection_time = 0;
            obj.performance_metrics.stability_score = 0;
        end
        
        function update_detection_statistics(obj, sgl_detected, flags)
            % 更新检测统计信息
            
            if sgl_detected
                obj.detection_stats.total_detections = ...
                    obj.detection_stats.total_detections + 1;
            end
            
            % 更新各准则成功率
            for i = 1:4
                if flags(i)
                    obj.detection_stats.criterion_success_rate(i) = ...
                        obj.detection_stats.criterion_success_rate(i) + 1;
                end
            end
            
            % 计算成功率（基于总检测次数）
            total_checks = length(obj.fisher_history);
            if total_checks > 0
                obj.detection_stats.criterion_success_rate = ...
                    obj.detection_stats.criterion_success_rate / total_checks;
            end
        end
        
        function log_detection_event(obj, fisher_info, innovation)
            % 记录检测事件（调试用）
            
            fprintf('[SGL检测] 时间: %.2f, 模式: SGL_DETECTED\n', obj.last_detection_time);
            fprintf('  Fisher信息迹: %.2e\n', trace(fisher_info));
            fprintf('  新息范数: %.4f\n', norm(innovation));
            fprintf('  检测标志: [%d, %d, %d, %d]\n', obj.detection_flags);
            
            % 计算理论饱和度
            saturation = obj.compute_saturation_ratio(fisher_info);
            fprintf('  理论饱和度: %.1f%%\n', saturation * 100);
        end
        
    end
    
    methods (Static)
        
        function run_tests()
            % 运行SGL检测器单元测试
            
            fprintf('=== SGL检测器单元测试 ===\n');
            
            % 创建测试配置
            config = struct();
            config.detection.tau1 = 0.1;
            config.detection.tau2 = 50;
            config.detection.tau3 = 0.05;
            config.detection.tau4 = 0.01;
            config.detection.window_length = 10;
            config.detection.confirm_count = 3;
            config.theory.C_sat = 2.5;
            config.theory.alpha = 0.67;
            config.theory.epsilon = 0.15;
            config.debug.verbose = false;
            
            % 测试1: 检测器初始化
            fprintf('测试1: 检测器初始化...');
            try
                detector = SGLDetector(config);
                fprintf(' 通过\n');
            catch ME
                fprintf(' 失败: %s\n', ME.message);
                return;
            end
            
            % 测试2: 理论界限计算
            fprintf('测试2: 理论界限计算...');
            bound_100 = detector.compute_theoretical_sgl_bound(100);
            bound_1000 = detector.compute_theoretical_sgl_bound(1000);
            if bound_1000 > bound_100 && bound_1000 < 10 * bound_100
                fprintf(' 通过 (n=100: %.3f, n=1000: %.3f)\n', bound_100, bound_1000);
            else
                fprintf(' 失败\n');
            end
            
            % 测试3: SGL检测逻辑
            fprintf('测试3: SGL检测逻辑...');
            
            % 模拟Fisher信息矩阵序列
            dim = 5;
            detection_results = [];
            
            for i = 1:20
                % 生成模拟Fisher信息矩阵
                if i <= 10
                    % 前期：信息快速增长
                    fisher = (i/10) * eye(dim) + 0.1 * randn(dim);
                    fisher = fisher' * fisher; % 确保正定
                else
                    % 后期：信息趋于饱和
                    fisher = eye(dim) + 0.01 * randn(dim);
                    fisher = fisher' * fisher;
                end
                
                innovation = 0.1 * randn(3, 1);
                S = eye(3);
                
                % 更新参数历史（模拟）
                if i == 1
                    detector.parameter_history = {};
                end
                theta = randn(dim, 1) * (1 - (i-1)/20); % 参数逐渐收敛
                detector.parameter_history{end+1} = theta;
                
                % 执行检测
                sgl_detected = detector.detect_sgl(fisher, innovation, S);
                detection_results = [detection_results, sgl_detected];
            end
            
            % 检查是否在后期检测到SGL
            late_detections = sum(detection_results(15:end));
            if late_detections > 0
                fprintf(' 通过 (后期检测次数: %d)\n', late_detections);
            else
                fprintf(' 失败 (未检测到SGL)\n');
            end
            
            % 测试4: 统计功能
            fprintf('测试4: 统计功能...');
            try
                stats = detector.get_detection_statistics();
                if isstruct(stats) && isfield(stats, 'total_detections')
                    fprintf(' 通过\n');
                else
                    fprintf(' 失败\n');
                end
            catch ME
                fprintf(' 失败: %s\n', ME.message);
            end
            
            fprintf('=== 测试完成 ===\n\n');
        end
        
    end
    
end

function str = logical_to_string(val)
% 将逻辑值转换为字符串
if val
    str = '满足';
else
    str = '不满足';
end
end