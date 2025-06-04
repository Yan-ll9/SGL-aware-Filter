function config = config_system(config_file)
% CONFIG_SYSTEM 系统配置管理
%
% 输入:
%   config_file - 配置文件路径（JSON格式）
%
% 输出:
%   config - 配置结构体
%
% 作者: Y
% 日期: 2025.6

if nargin < 1 || isempty(config_file)
    config_file = 'config/filter_params.json';
end

% 检查配置文件是否存在
if ~exist(config_file, 'file')
    fprintf('配置文件不存在: %s\n', config_file);
    fprintf('使用默认配置...\n');
    config = get_default_config();
    return;
end

try
    % 读取JSON配置文件
    config_text = fileread(config_file);
    config = jsondecode(config_text);
    
    % 验证和补充配置
    config = validate_and_complete_config(config);
    
    fprintf('配置文件加载成功: %s\n', config_file);
    
catch ME
    fprintf('配置文件解析失败: %s\n', ME.message);
    fprintf('使用默认配置...\n');
    config = get_default_config();
end

end

function config = get_default_config()
% 获取默认配置

config = struct();

%% 系统基本参数
config.system.name = 'SGL-aware GNSS/INS Filter';
config.system.version = '1.0.0';
config.system.author = 'Research Team';

%% 数据配置
config.data.imu_freq = 100;           % IMU频率 (Hz)
config.data.gnss_freq = 10;           % GNSS频率 (Hz)
config.data.sync.tolerance = 0.01;    % 同步容差 (s)
config.data.sync.method = 'nearest';  % 同步方法

% 数据格式配置
config.data.formats.imu = 'time,ax,ay,az,gx,gy,gz';
config.data.formats.gnss = 'time,x,y,z,vx,vy,vz,std_pos,std_vel';
config.data.formats.reference = 'time,x,y,z,vx,vy,vz,roll,pitch,yaw';

%% 滤波器基本配置
config.filter.type = 'SGL-aware';
config.filter.state_dim = 17;         % 状态维度

% 状态定义: [pos(3), vel(3), att_log(3), ba(3), bg(3), clk(2)]
config.filter.state_names = {
    'pos_x', 'pos_y', 'pos_z', ...
    'vel_x', 'vel_y', 'vel_z', ...
    'att_x', 'att_y', 'att_z', ...
    'ba_x', 'ba_y', 'ba_z', ...
    'bg_x', 'bg_y', 'bg_z', ...
    'clk_bias', 'clk_drift'
};

% 观测维度配置
config.filter.obs_dim_gnss = 6;       % GNSS观测维度 (位置+速度)

%% SGL理论参数
config.sgl.enable = true;
config.sgl.detection.tau1 = 0.1;      % 谱饱和准则阈值
config.sgl.detection.tau2 = 50;       % 条件数稳定准则阈值  
config.sgl.detection.tau3 = 0.05;     % 信息增量衰减准则阈值
config.sgl.detection.tau4 = 0.01;     % 参数收敛准则阈值
config.sgl.detection.window_length = 100;  % 观察窗口长度
config.sgl.detection.confirm_count = 5;    % 确认阈值

% SGL理论常数（根据文章理论推导）
config.sgl.theory.C_sat = 2.5;        % 饱和常数
config.sgl.theory.alpha = 0.67;       % 几何约束参数
config.sgl.theory.epsilon = 0.15;     % 几何约束参数

%% 四模式控制参数
config.modes.active.alpha_base = 0.01;     % ACTIVE模式基础步长
config.modes.observe.alpha_base = 0.003;   % OBSERVE模式基础步长
config.modes.frozen.alpha_base = 0;        % FROZEN模式步长
config.modes.reactivate.alpha_base = 0.02; % REACTIVATE模式步长

config.modes.switch.gamma_min = 0.1;       % 最小残差阈值
config.modes.switch.chi2_threshold = 15.5; % 卡方检验阈值(95%置信度,6自由度)

%% 初始状态参数
config.init.pos_var = 100;            % 初始位置方差 (m^2)
config.init.vel_var = 10;             % 初始速度方差 (m^2/s^2)
config.init.att_var = (5*pi/180)^2;   % 初始姿态方差 (rad^2)
config.init.ba_var = (0.1)^2;         % 初始加速度计偏差方差 (m^2/s^4)
config.init.bg_var = (1*pi/180)^2;    % 初始陀螺仪偏差方差 (rad^2/s^2)
config.init.clock_var = 100;          % 初始时钟方差 (m^2)

%% 过程噪声参数（初始值，会被自适应调整）
config.noise.process.pos_noise = 0;              % 位置过程噪声
config.noise.process.vel_noise = 0.1;            % 速度过程噪声 (m/s^2)
config.noise.process.att_noise = 1*pi/180;       % 姿态过程噪声 (rad/s)
config.noise.process.ba_noise = 0.01;            % 加速度计偏差噪声 (m/s^3)
config.noise.process.bg_noise = 0.1*pi/180;      % 陀螺仪偏差噪声 (rad/s^2)
config.noise.process.clock_bias_noise = 1;       % 时钟偏差噪声 (m/s)
config.noise.process.clock_drift_noise = 0.1;    % 时钟漂移噪声 (m/s^2)

% IMU相关时间常数
config.noise.correlation_time.acc_bias = 3600;   % 加速度计偏差相关时间 (s)
config.noise.correlation_time.gyro_bias = 3600;  % 陀螺仪偏差相关时间 (s)

%% 观测噪声参数（初始值，会被自适应调整）
config.noise.observation.gnss_pos = 5;           % GNSS位置观测噪声 (m)
config.noise.observation.gnss_vel = 0.1;         % GNSS速度观测噪声 (m/s)

%% 自适应参数更新配置
config.adaptive.enable = true;
config.adaptive.update_interval = 10;            % 更新间隔（历元数）
config.adaptive.forgetting_factor = 0.95;       % 遗忘因子
config.adaptive.min_samples = 50;               % 最小样本数

% 自然梯度参数
config.adaptive.natural_gradient.C_step = 1.0;  % 步长控制常数

%% 几何参数
config.geometry.gravity = [0; 0; -9.8061];       % 重力加速度 (m/s^2)
config.geometry.earth_rotation = 7.2921159e-5;   % 地球自转角速度 (rad/s)
config.geometry.enable_earth_rotation = false;   % 是否考虑地球自转

%% 数值积分参数
config.integration.method = 'rk4';               % 积分方法
config.integration.max_step = 0.01;              % 最大积分步长 (s)
config.integration.tolerance = 1e-6;             % 积分容差

%% 异常检测参数
config.outlier_detection.enable = true;
config.outlier_detection.method = 'chi2';        % 异常检测方法
config.outlier_detection.threshold = 13.8;      % 卡方阈值(99%置信度,6自由度)
config.outlier_detection.max_outliers = 0.1;    % 最大异常比例

%% 数据质量控制
config.quality_control.min_satellites = 4;      % 最少卫星数
config.quality_control.max_hdop = 10;           % 最大HDOP
config.quality_control.min_elevation = 10;      % 最小高度角 (deg)
config.quality_control.max_age = 30;            % 最大数据龄期 (s)

%% 可视化配置
config.visualization.enable = true;
config.visualization.update_interval = 100;     % 可视化更新间隔
config.visualization.save_plots = true;
config.visualization.plot_types = {
    'trajectory', 'errors', 'sgl_analysis', 'mode_switching', ...
    'noise_adaptation', 'fisher_info', 'covariance'
};

%% 后处理配置
config.post_process.enable_smoothing = false;
config.post_process.smoothing.method = 'rts';   % RTS平滑
config.post_process.smoothing.lag = 10;         % 平滑滞后

%% 输出配置
config.output.save_full_states = true;
config.output.save_covariances = true;
config.output.save_innovations = true;
config.output.save_fisher_info = false;         % Fisher信息矩阵较大，可选保存
config.output.decimation_factor = 1;            % 输出抽取因子

%% 调试配置
config.debug.verbose = false;
config.debug.save_intermediate = false;
config.debug.check_constraints = true;          % 检查流形约束
config.debug.timing = false;                    % 性能计时

%% 并行计算配置
config.parallel.enable = false;
config.parallel.num_workers = 4;
config.parallel.use_gpu = false;

end

function config = validate_and_complete_config(config)
% 验证并补全配置

% 获取默认配置作为参考
default_config = get_default_config();

% 递归补全缺失字段
config = complete_config_fields(config, default_config);

% 验证关键参数
config = validate_config_parameters(config);

% 计算派生参数
config = compute_derived_parameters(config);

end

function config = complete_config_fields(config, default_config)
% 递归补全配置字段

field_names = fieldnames(default_config);

for i = 1:length(field_names)
    field_name = field_names{i};
    
    if ~isfield(config, field_name)
        % 字段不存在，使用默认值
        config.(field_name) = default_config.(field_name);
    elseif isstruct(default_config.(field_name))
        % 递归处理结构体字段
        if ~isstruct(config.(field_name))
            config.(field_name) = default_config.(field_name);
        else
            config.(field_name) = complete_config_fields(config.(field_name), default_config.(field_name));
        end
    end
end

end

function config = validate_config_parameters(config)
% 验证配置参数的合理性

% 验证频率参数
if config.data.imu_freq <= 0
    warning('IMU频率无效，使用默认值100Hz');
    config.data.imu_freq = 100;
end

if config.data.gnss_freq <= 0
    warning('GNSS频率无效，使用默认值10Hz');
    config.data.gnss_freq = 10;
end

% 验证SGL参数
if config.sgl.detection.tau1 <= 0 || config.sgl.detection.tau1 >= 1
    warning('SGL tau1参数无效，使用默认值0.1');
    config.sgl.detection.tau1 = 0.1;
end

if config.sgl.theory.alpha <= 0 || config.sgl.theory.alpha >= 1
    warning('SGL alpha参数无效，使用默认值0.67');
    config.sgl.theory.alpha = 0.67;
end

% 验证噪声参数
noise_fields = fieldnames(config.noise.process);
for i = 1:length(noise_fields)
    field_name = noise_fields{i};
    if config.noise.process.(field_name) < 0
        warning('过程噪声参数 %s 为负值，设为默认值', field_name);
        default_config = get_default_config();
        config.noise.process.(field_name) = default_config.noise.process.(field_name);
    end
end

% 验证协方差参数
if any(diag([config.init.pos_var, config.init.vel_var, config.init.att_var, ...
             config.init.ba_var, config.init.bg_var, config.init.clock_var]) <= 0)
    warning('初始协方差参数存在非正值');
end

end

function config = compute_derived_parameters(config)
% 计算派生参数

% 计算时间间隔
config.dt.imu = 1 / config.data.imu_freq;
config.dt.gnss = 1 / config.data.gnss_freq;

% 计算离散化过程噪声协方差矩阵
config.computed.Q_discrete = compute_discrete_process_noise(config);

% 计算观测噪声协方差矩阵
config.computed.R_gnss = diag([
    config.noise.observation.gnss_pos^2 * ones(3, 1);
    config.noise.observation.gnss_vel^2 * ones(3, 1)
]);

% 计算初始协方差矩阵
config.computed.P0 = blkdiag(
    config.init.pos_var * eye(3), ...
    config.init.vel_var * eye(3), ...
    config.init.att_var * eye(3), ...
    config.init.ba_var * eye(3), ...
    config.init.bg_var * eye(3), ...
    config.init.clock_var * eye(2)
);

% SGL理论相关计算
config.computed.sgl_saturation_bound = config.sgl.theory.C_sat / config.sgl.theory.epsilon;

% 计算几何常数
config.computed.gravity_norm = norm(config.geometry.gravity);

end

function Q_discrete = compute_discrete_process_noise(config)
% 计算离散化过程噪声协方差矩阵

dt = config.dt.imu;

% 连续时间过程噪声功率谱密度
S_pos = config.noise.process.pos_noise^2;
S_vel = config.noise.process.vel_noise^2;
S_att = config.noise.process.att_noise^2;
S_ba = config.noise.process.ba_noise^2;
S_bg = config.noise.process.bg_noise^2;
S_clk_bias = config.noise.process.clock_bias_noise^2;
S_clk_drift = config.noise.process.clock_drift_noise^2;

% 离散化（Van Loan方法）
% 位置（积分白噪声模型）
Q_pos = S_pos * dt * eye(3);

% 速度（白噪声模型）
Q_vel = S_vel * dt * eye(3);

% 姿态（随机游走模型）
Q_att = S_att * dt * eye(3);

% IMU偏差（一阶马尔可夫过程）
tau_a = config.noise.correlation_time.acc_bias;
tau_g = config.noise.correlation_time.gyro_bias;

if tau_a > 0
    Q_ba = S_ba * (1 - exp(-2*dt/tau_a)) / 2 * eye(3);
else
    Q_ba = S_ba * dt * eye(3);
end

if tau_g > 0
    Q_bg = S_bg * (1 - exp(-2*dt/tau_g)) / 2 * eye(3);
else
    Q_bg = S_bg * dt * eye(3);
end

% 时钟参数（二阶随机游走模型）
Q_clk = [dt^3/3, dt^2/2; dt^2/2, dt] * S_clk_drift + ...
        [dt, 0; 0, 0] * S_clk_bias;

% 组装完整的过程噪声协方差矩阵
Q_discrete = blkdiag(Q_pos, Q_vel, Q_att, Q_ba, Q_bg, Q_clk);

end

function save_config(config, filename)
% 保存配置到文件

if nargin < 2
    filename = sprintf('config/saved_config_%s.json', datestr(now, 'yyyymmdd_HHMMSS'));
end

% 确保目录存在
[filepath, ~, ~] = fileparts(filename);
if ~exist(filepath, 'dir')
    mkdir(filepath);
end

try
    % 转换为JSON并保存
    config_json = jsonencode(config, 'PrettyPrint', true);
    
    fid = fopen(filename, 'w');
    if fid == -1
        error('无法创建配置文件: %s', filename);
    end
    
    fprintf(fid, '%s', config_json);
    fclose(fid);
    
    fprintf('配置已保存到: %s\n', filename);
    
catch ME
    warning('配置保存失败: %s', ME.message);
end

end

function display_config_summary(config)
% 显示配置摘要

fprintf('\n=== SGL-aware滤波器配置摘要 ===\n');
fprintf('系统: %s v%s\n', config.system.name, config.system.version);
fprintf('作者: %s\n', config.system.author);

fprintf('\n--- 数据配置 ---\n');
fprintf('IMU频率: %d Hz\n', config.data.imu_freq);
fprintf('GNSS频率: %d Hz\n', config.data.gnss_freq);
fprintf('同步容差: %.3f s\n', config.data.sync.tolerance);

fprintf('\n--- 滤波器配置 ---\n');
fprintf('类型: %s\n', config.filter.type);
fprintf('状态维度: %d\n', config.filter.state_dim);
fprintf('SGL使能: %s\n', logical_to_string(config.sgl.enable));

fprintf('\n--- SGL参数 ---\n');
fprintf('检测阈值: [%.3f, %.1f, %.3f, %.3f]\n', ...
        config.sgl.detection.tau1, config.sgl.detection.tau2, ...
        config.sgl.detection.tau3, config.sgl.detection.tau4);
fprintf('理论参数: C_sat=%.2f, α=%.2f, ε=%.2f\n', ...
        config.sgl.theory.C_sat, config.sgl.theory.alpha, config.sgl.theory.epsilon);

fprintf('\n--- 自适应配置 ---\n');
fprintf('自适应使能: %s\n', logical_to_string(config.adaptive.enable));
fprintf('更新间隔: %d 历元\n', config.adaptive.update_interval);
fprintf('遗忘因子: %.3f\n', config.adaptive.forgetting_factor);

fprintf('\n--- 质量控制 ---\n');
fprintf('最少卫星数: %d\n', config.quality_control.min_satellites);
fprintf('最大HDOP: %.1f\n', config.quality_control.max_hdop);
fprintf('异常检测: %s\n', logical_to_string(config.outlier_detection.enable));

fprintf('===============================\n\n');

end

function str = logical_to_string(val)
% 将逻辑值转换为字符串
if val
    str = '是';
else
    str = '否';
end
end

% 如果直接运行此文件，显示默认配置
if ~nargout
    default_config = get_default_config();
    display_config_summary(default_config);
    save_config(default_config, 'config/default_config.json');
end