function main_sgl_filter(data_path, config_file, varargin)
% MAIN_SGL_FILTER SGL-aware自适应滤波主程序
%
% 输入:
%   data_path   - 数据文件路径或实时数据源配置
%   config_file - 配置文件路径 (可选)
%   varargin    - 可选参数 ('realtime', 'plot', 'save_results', etc.)
%
% 示例:
%   main_sgl_filter('data_samples/real_data/urban_canyon/', 'config/default.json')
%   main_sgl_filter('data_samples/real_data/urban_canyon/', [], 'plot', 'save_results')
%
% 作者: Y
% 日期: 2025.6

%% 初始化
addpath(genpath(pwd)); % 添加所有子目录到路径
rng(42); % 设置随机种子，确保结果可重复

% 解析输入参数
p = inputParser;
addRequired(p, 'data_path', @(x) ischar(x) || isstring(x));
addOptional(p, 'config_file', 'config/filter_params.json', @(x) ischar(x) || isstring(x));
addParameter(p, 'realtime', false, @islogical);
addParameter(p, 'plot', true, @islogical);
addParameter(p, 'save_results', true, @islogical);
addParameter(p, 'verbose', true, @islogical);
addParameter(p, 'output_dir', 'results/', @(x) ischar(x) || isstring(x));

parse(p, data_path, config_file, varargin{:});
args = p.Results;

if args.verbose
    fprintf('=== SGL-aware自适应滤波系统启动 ===\n');
    fprintf('数据路径: %s\n', args.data_path);
    fprintf('配置文件: %s\n', args.config_file);
    fprintf('实时模式: %s\n', string(args.realtime));
end

%% 加载配置
try
    config = config_system(args.config_file);
    if args.verbose
        fprintf('配置加载成功\n');
    end
catch ME
    error('配置加载失败: %s', ME.message);
end

%% 初始化数据接口
try
    if args.realtime
        % 实时数据接口
        data_interface = RealTimeDataInterface(args.data_path, config.data);
        if args.verbose
            fprintf('实时数据接口初始化成功\n');
        end
    else
        % 离线数据加载
        data_loader = DataLoader(args.data_path, config.data);
        [imu_data, gnss_data, ref_data] = data_loader.load_all();
        if args.verbose
            fprintf('数据加载成功: IMU数据点 %d, GNSS数据点 %d\n', ...
                    size(imu_data, 1), size(gnss_data, 1));
        end
    end
catch ME
    error('数据加载失败: %s', ME.message);
end

%% 初始化SGL-aware滤波器
try
    sgl_filter = SGLAwareFilter(config);
    if args.verbose
        fprintf('SGL-aware滤波器初始化成功\n');
    end
catch ME
    error('滤波器初始化失败: %s', ME.message);
end

%% 初始化可视化和结果保存
if args.plot
    visualizer = create_visualizer(config.visualization);
end

if args.save_results
    if ~exist(args.output_dir, 'dir')
        mkdir(args.output_dir);
    end
    result_saver = ResultSaver(args.output_dir);
end

%% 主处理循环
if args.verbose
    fprintf('\n=== 开始滤波处理 ===\n');
    tic;
end

try
    if args.realtime
        % 实时处理模式
        results = process_realtime(sgl_filter, data_interface, config, args);
    else
        % 离线处理模式
        results = process_offline(sgl_filter, imu_data, gnss_data, ref_data, config, args);
    end
    
    processing_time = toc;
    if args.verbose
        fprintf('滤波处理完成，耗时 %.2f 秒\n', processing_time);
    end
    
catch ME
    fprintf('处理过程中发生错误: %s\n', ME.message);
    fprintf('错误位置: %s (第 %d 行)\n', ME.stack(1).file, ME.stack(1).line);
    rethrow(ME);
end

%% 结果分析和可视化
if args.verbose
    fprintf('\n=== 结果分析 ===\n');
end

try
    % 性能统计
    perf_stats = analyze_performance(results, config);
    
    if args.verbose
        display_performance_summary(perf_stats);
    end
    
    % 可视化
    if args.plot
        create_plots(results, perf_stats, visualizer, config);
    end
    
    % 保存结果
    if args.save_results
        save_all_results(results, perf_stats, result_saver, config);
    end
    
catch ME
    warning('结果分析过程中发生错误: %s', ME.message);
end

if args.verbose
    fprintf('\n=== 处理完成 ===\n');
end

end

%% 离线处理函数
function results = process_offline(sgl_filter, imu_data, gnss_data, ref_data, config, args)
    
    % 数据同步和预处理
    data_sync = DataSynchronizer(config.data.sync);
    [imu_sync, gnss_sync, ref_sync] = data_sync.synchronize(imu_data, gnss_data, ref_data);
    
    % 初始化结果存储
    num_epochs = size(imu_sync, 1);
    results = initialize_results_structure(num_epochs, config);
    
    % 状态初始化
    [X0, P0] = initialize_state(imu_sync(1:100, :), gnss_sync(1, :), config);
    sgl_filter.initialize(X0, P0);
    
    % 主处理循环
    gnss_idx = 1;
    update_interval = max(1, floor(num_epochs / 100)); % 更新进度的间隔
    
    for k = 1:num_epochs
        % 显示进度
        if args.verbose && mod(k, update_interval) == 0
            fprintf('处理进度: %d/%d (%.1f%%)\n', k, num_epochs, 100*k/num_epochs);
        end
        
        % IMU预测
        imu_meas = imu_sync(k, :);
        sgl_filter.predict(imu_meas);
        
        % GNSS更新（如果有新观测）
        if gnss_idx <= size(gnss_sync, 1) && ...
           abs(gnss_sync(gnss_idx, 1) - imu_sync(k, 1)) < config.data.sync.tolerance
            
            gnss_meas = gnss_sync(gnss_idx, :);
            sgl_filter.update(gnss_meas);
            gnss_idx = gnss_idx + 1;
        end
        
        % 保存结果
        results = save_epoch_results(results, k, sgl_filter, ref_sync, k);
    end
    
    % 后处理和平滑
    results = post_process_results(results, config);
end

%% 实时处理函数
function results = process_realtime(sgl_filter, data_interface, config, args)
    
    results = [];
    k = 1;
    
    % 初始化
    data_interface.start();
    
    try
        while data_interface.has_data()
            % 获取新数据
            [data_type, measurement, timestamp] = data_interface.get_next();
            
            switch data_type
                case 'IMU'
                    sgl_filter.predict(measurement);
                case 'GNSS'
                    sgl_filter.update(measurement);
            end
            
            % 定期保存结果
            if mod(k, 10) == 0
                current_state = sgl_filter.get_current_state();
                results(end+1, :) = [timestamp, current_state]; %#ok<AGROW>
            end
            
            k = k + 1;
        end
    catch ME
        data_interface.stop();
        rethrow(ME);
    end
    
    data_interface.stop();
end

%% 辅助函数
function results = initialize_results_structure(num_epochs, config)
    % 初始化结果数据结构
    state_dim = config.filter.state_dim;
    
    results.timestamps = zeros(num_epochs, 1);
    results.states = zeros(num_epochs, state_dim);
    results.covariances = zeros(state_dim, state_dim, num_epochs);
    results.innovations = cell(num_epochs, 1);
    results.sgl_modes = zeros(num_epochs, 1);
    results.fisher_info = cell(num_epochs, 1);
    results.noise_params = cell(num_epochs, 1);
    results.performance_metrics = zeros(num_epochs, 10); % 预留10个性能指标
end

function [X0, P0] = initialize_state(imu_init, gnss_init, config)
    % 状态初始化
    
    % 位置初始化（从GNSS）
    pos_init = gnss_init(2:4)'; % 假设格式: [time, x, y, z, ...]
    
    % 速度初始化（从GNSS或设为零）
    if size(gnss_init, 2) >= 7
        vel_init = gnss_init(5:7)';
    else
        vel_init = zeros(3, 1);
    end
    
    % 姿态初始化（从IMU静态对准）
    acc_init = mean(imu_init(:, 2:4), 1)'; % 假设格式: [time, ax, ay, az, gx, gy, gz]
    gyro_init = mean(imu_init(:, 5:7), 1)';
    
    % 重力对准获得初始姿态
    R_init = gravity_alignment(acc_init, config.gravity);
    
    % IMU偏差初始化
    ba_init = acc_init + R_init' * config.gravity; % 加速度计偏差
    bg_init = gyro_init; % 陀螺仪偏差
    
    % 时钟偏差初始化
    clock_init = [0; 0]; % [偏差; 漂移]
    
    % 组装状态向量
    X0 = [pos_init; vel_init; R_init(:); ba_init; bg_init; clock_init];
    
    % 初始协方差
    P0 = blkdiag(
        config.init.pos_var * eye(3), ...      % 位置
        config.init.vel_var * eye(3), ...      % 速度  
        config.init.att_var * eye(3), ...      % 姿态（李代数表示）
        config.init.ba_var * eye(3), ...       % 加速度计偏差
        config.init.bg_var * eye(3), ...       % 陀螺仪偏差
        config.init.clock_var * eye(2) ...     % 时钟参数
    );
end

function R = gravity_alignment(acc_meas, gravity_ref)
    % 基于重力矢量的姿态初始对准
    
    % 归一化
    acc_norm = acc_meas / norm(acc_meas);
    grav_norm = gravity_ref / norm(gravity_ref);
    
    % 计算旋转轴和角度
    v = cross(acc_norm, -grav_norm); % 注意加速度计测量的是-g
    s = norm(v);
    c = dot(acc_norm, -grav_norm);
    
    if s < 1e-6
        R = eye(3);
    else
        vx = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
        R = eye(3) + vx + vx^2 * (1-c)/s^2;
    end
end

function perf_stats = analyze_performance(results, config)
    % 性能分析
    
    perf_stats = struct();
    
    % 计算位置误差
    if isfield(results, 'reference')
        pos_errors = results.states(:, 1:3) - results.reference(:, 1:3);
        perf_stats.pos_rmse = sqrt(mean(sum(pos_errors.^2, 2)));
        perf_stats.pos_std = std(sqrt(sum(pos_errors.^2, 2)));
        
        % 计算速度误差
        vel_errors = results.states(:, 4:6) - results.reference(:, 4:6);
        perf_stats.vel_rmse = sqrt(mean(sum(vel_errors.^2, 2)));
        
        % 计算姿态误差
        att_errors = compute_attitude_errors(results.states(:, 7:15), results.reference(:, 7:15));
        perf_stats.att_rmse = sqrt(mean(sum(att_errors.^2, 2)));
    end
    
    % SGL统计
    mode_counts = histcounts(results.sgl_modes, 1:5);
    perf_stats.mode_distribution = mode_counts / sum(mode_counts);
    
    % 计算效率指标
    perf_stats.avg_compute_load = mean(results.performance_metrics(:, 1));
    perf_stats.convergence_time = find_convergence_time(results);
end

function display_performance_summary(perf_stats)
    % 显示性能摘要
    
    fprintf('=== 性能统计摘要 ===\n');
    if isfield(perf_stats, 'pos_rmse')
        fprintf('位置RMSE: %.3f m\n', perf_stats.pos_rmse);
        fprintf('速度RMSE: %.3f m/s\n', perf_stats.vel_rmse);
        fprintf('姿态RMSE: %.3f deg\n', rad2deg(perf_stats.att_rmse));
    end
    
    fprintf('模式分布:\n');
    mode_names = {'ACTIVE', 'OBSERVE', 'FROZEN', 'REACTIVATE'};
    for i = 1:4
        fprintf('  %s: %.1f%%\n', mode_names{i}, perf_stats.mode_distribution(i)*100);
    end
    
    fprintf('平均计算负荷: %.1f%%\n', perf_stats.avg_compute_load);
    fprintf('收敛时间: %.1f s\n', perf_stats.convergence_time);
end

function create_plots(results, perf_stats, visualizer, config)
    % 创建可视化图表
    
    % 轨迹图
    visualizer.plot_trajectory(results.states(:, 1:3), results.reference(:, 1:3));
    
    % 误差时间序列
    visualizer.plot_errors(results);
    
    % SGL分析图
    visualizer.plot_sgl_analysis(results);
    
    % 模式切换图
    visualizer.plot_mode_switching(results.sgl_modes, results.timestamps);
end

function visualizer = create_visualizer(viz_config)
    % 创建可视化器
    visualizer = struct();
    visualizer.config = viz_config;
    
    % 添加绘图方法
    visualizer.plot_trajectory = @plot_trajectory_impl;
    visualizer.plot_errors = @plot_errors_impl;
    visualizer.plot_sgl_analysis = @plot_sgl_analysis_impl;
    visualizer.plot_mode_switching = @plot_mode_switching_impl;
end

function plot_trajectory_impl(estimated_pos, reference_pos)
    figure('Name', '轨迹对比', 'Position', [100, 100, 800, 600]);
    plot3(estimated_pos(:,1), estimated_pos(:,2), estimated_pos(:,3), 'b-', 'LineWidth', 2);
    hold on;
    if ~isempty(reference_pos)
        plot3(reference_pos(:,1), reference_pos(:,2), reference_pos(:,3), 'r--', 'LineWidth', 1);
        legend('SGL-aware估计', '参考轨迹', 'Location', 'best');
    end
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    title('三维轨迹对比');
    grid on; axis equal;
end

function plot_errors_impl(results)
    if ~isfield(results, 'reference')
        return;
    end
    
    figure('Name', '误差分析', 'Position', [200, 200, 1200, 800]);
    
    % 位置误差
    subplot(3,1,1);
    pos_errors = results.states(:, 1:3) - results.reference(:, 1:3);
    plot(results.timestamps, pos_errors);
    ylabel('位置误差 (m)');
    legend('X', 'Y', 'Z');
    title('位置误差时间序列');
    grid on;
    
    % 速度误差
    subplot(3,1,2);
    vel_errors = results.states(:, 4:6) - results.reference(:, 4:6);
    plot(results.timestamps, vel_errors);
    ylabel('速度误差 (m/s)');
    legend('Vx', 'Vy', 'Vz');
    title('速度误差时间序列');
    grid on;
    
    % 姿态误差
    subplot(3,1,3);
    att_errors = compute_attitude_errors(results.states(:, 7:15), results.reference(:, 7:15));
    plot(results.timestamps, rad2deg(att_errors));
    ylabel('姿态误差 (deg)');
    xlabel('时间 (s)');
    legend('Roll', 'Pitch', 'Yaw');
    title('姿态误差时间序列');
    grid on;
end

function plot_sgl_analysis_impl(results)
    figure('Name', 'SGL分析', 'Position', [300, 300, 1200, 600]);
    
    % Fisher信息演化
    subplot(2,2,1);
    fisher_traces = cellfun(@(x) trace(x), results.fisher_info);
    plot(results.timestamps, fisher_traces);
    ylabel('Fisher信息迹');
    title('Fisher信息积累');
    grid on;
    
    % SGL检测指标
    subplot(2,2,2);
    % 这里需要从results中提取SGL检测指标
    % plot(results.timestamps, sgl_criteria);
    title('SGL检测准则');
    grid on;
    
    % 噪声参数演化
    subplot(2,2,3);
    % 绘制关键噪声参数的变化
    title('噪声参数自适应');
    grid on;
    
    % 计算负荷
    subplot(2,2,4);
    plot(results.timestamps, results.performance_metrics(:, 1));
    ylabel('计算负荷 (%)');
    xlabel('时间 (s)');
    title('计算负荷变化');
    grid on;
end

function plot_mode_switching_impl(modes, timestamps)
    figure('Name', '模式切换', 'Position', [400, 400, 1000, 300]);
    
    % 创建模式切换图
    stairs(timestamps, modes, 'LineWidth', 2);
    ylim([0.5, 4.5]);
    yticks(1:4);
    yticklabels({'ACTIVE', 'OBSERVE', 'FROZEN', 'REACTIVATE'});
    xlabel('时间 (s)');
    ylabel('运行模式');
    title('SGL-aware四模式切换时序');
    grid on;
end

function att_errors = compute_attitude_errors(est_rotmats, ref_rotmats)
    % 计算姿态误差（欧拉角形式）
    
    n = size(est_rotmats, 1);
    att_errors = zeros(n, 3);
    
    for i = 1:n
        R_est = reshape(est_rotmats(i, :), 3, 3);
        R_ref = reshape(ref_rotmats(i, :), 3, 3);
        
        % 相对旋转矩阵
        R_err = R_ref' * R_est;
        
        % 转换为欧拉角误差
        att_errors(i, :) = rotm2eul(R_err, 'ZYX');
    end
end

function conv_time = find_convergence_time(results)
    % 寻找收敛时间
    
    % 简单方法：寻找模式首次切换到FROZEN的时间
    frozen_idx = find(results.sgl_modes == 3, 1, 'first');
    if ~isempty(frozen_idx)
        conv_time = results.timestamps(frozen_idx);
    else
        conv_time = inf;
    end
end

function results = save_epoch_results(results, k, sgl_filter, ref_data, epoch_idx)
    % 保存单个历元的结果
    
    current_state = sgl_filter.get_current_state();
    current_cov = sgl_filter.get_current_covariance();
    
    results.timestamps(k) = current_state.timestamp;
    results.states(k, :) = current_state.X;
    results.covariances(:, :, k) = current_cov;
    results.sgl_modes(k) = sgl_filter.get_current_mode();
    results.fisher_info{k} = sgl_filter.get_fisher_info();
    results.noise_params{k} = sgl_filter.get_noise_params();
    
    % 保存参考数据（如果有）
    if ~isempty(ref_data) && epoch_idx <= size(ref_data, 1)
        results.reference(k, :) = ref_data(epoch_idx, 2:end); % 假设第一列是时间
    end
    
    % 保存性能指标
    perf_metrics = sgl_filter.get_performance_metrics();
    results.performance_metrics(k, :) = perf_metrics;
end

function results = post_process_results(results, config)
    % 后处理结果
    
    % 平滑处理（可选）
    if config.post_process.enable_smoothing
        results = apply_smoothing(results, config.post_process.smoothing);
    end
    
    % 计算额外的统计量
    results = compute_additional_statistics(results);
end

function save_all_results(results, perf_stats, result_saver, config)
    % 保存所有结果
    
    timestamp_str = datestr(now, 'yyyymmdd_HHMMSS');
    
    % 保存主要结果
    save([result_saver.output_dir, 'results_', timestamp_str, '.mat'], 'results', 'perf_stats', 'config');
    
    % 保存CSV格式的轨迹
    trajectory_table = table(results.timestamps, results.states(:,1), results.states(:,2), results.states(:,3), ...
                           'VariableNames', {'Time', 'X', 'Y', 'Z'});
    writetable(trajectory_table, [result_saver.output_dir, 'trajectory_', timestamp_str, '.csv']);
    
    % 保存性能报告
    save_performance_report(perf_stats, [result_saver.output_dir, 'performance_', timestamp_str, '.txt']);
end

function save_performance_report(perf_stats, filename)
    % 保存性能报告
    
    fid = fopen(filename, 'w');
    fprintf(fid, 'SGL-aware自适应滤波性能报告\n');
    fprintf(fid, '生成时间: %s\n\n', datestr(now));
    
    if isfield(perf_stats, 'pos_rmse')
        fprintf(fid, '=== 精度指标 ===\n');
        fprintf(fid, '位置RMSE: %.3f m\n', perf_stats.pos_rmse);
        fprintf(fid, '速度RMSE: %.3f m/s\n', perf_stats.vel_rmse);
        fprintf(fid, '姿态RMSE: %.3f deg\n', rad2deg(perf_stats.att_rmse));
    end
    
    fprintf(fid, '\n=== 效率指标 ===\n');
    fprintf(fid, '平均计算负荷: %.1f%%\n', perf_stats.avg_compute_load);
    fprintf(fid, '收敛时间: %.1f s\n', perf_stats.convergence_time);
    
    fprintf(fid, '\n=== 模式分布 ===\n');
    mode_names = {'ACTIVE', 'OBSERVE', 'FROZEN', 'REACTIVATE'};
    for i = 1:4
        fprintf(fid, '%s: %.1f%%\n', mode_names{i}, perf_stats.mode_distribution(i)*100);
    end
    
    fclose(fid);
end

function result_saver = ResultSaver(output_dir)
    % 结果保存器构造函数
    result_saver = struct();
    result_saver.output_dir = output_dir;
end