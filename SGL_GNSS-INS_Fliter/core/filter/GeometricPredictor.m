classdef GeometricPredictor < handle
% GEOMETRICPREDICTOR 几何预测器类
%
% 
% 核心功能：
% - 复合流形上的Munthe-Kaas数值积分
% - 几何约束保持的状态传播
% - 协方差预测的几何线性化
% - IMU动力学建模
%
% 
%
% 作者: Y
% 日期: 2025.6

    properties (Access = private)
        config          % 配置参数
        gravity         % 重力向量 [m/s²]
        earth_rotation  % 地球自转角速度 [rad/s]
        
        % 积分参数
        integration_method  % 积分方法
        max_step           % 最大积分步长
        tolerance          % 积分容差
        
        % 性能统计
        prediction_times   % 预测计算时间
        constraint_violations  % 约束违反次数
        
        % 调试信息
        debug_enabled     % 是否启用调试
        last_prediction   % 上次预测结果
    end
    
    methods (Access = public)
        
        function obj = GeometricPredictor(config)
            % 构造函数
            %
            % 输入:
            %   config - 配置结构体
            
            obj.config = config;
            
            % 初始化物理参数
            obj.gravity = config.geometry.gravity;
            obj.earth_rotation = config.geometry.earth_rotation;
            
            % 初始化积分参数
            obj.integration_method = config.integration.method;
            obj.max_step = config.integration.max_step;
            obj.tolerance = config.integration.tolerance;
            
            % 初始化性能统计
            obj.prediction_times = [];
            obj.constraint_violations = 0;
            
            % 调试设置
            obj.debug_enabled = config.debug.verbose;
        end
        
        function [X_pred, P_pred] = predict(obj, X_current, P_current, acc_meas, gyro_meas, dt, Q)
            % 主预测函数
            %
            % 输入:
            %   X_current - 17x1 当前状态向量
            %   P_current - 15x15 当前协方差矩阵(切空间)
            %   acc_meas - 3x1 加速度计测量 [m/s²]
            %   gyro_meas - 3x1 陀螺仪测量 [rad/s]
            %   dt - 时间间隔 [s]
            %   Q - 过程噪声协方差矩阵
            %
            % 输出:
            %   X_pred - 17x1 预测状态向量
            %   P_pred - 15x15 预测协方差矩阵
            
            tic; % 开始计时
            
            % 输入验证
            if ~ManifoldOperations.is_valid_state(X_current)
                warning('当前状态无效，正在修正...');
                X_current = ManifoldOperations.enforce_constraints(X_current);
            end
            
            % 状态预测 - 使用几何数值积分
            X_pred = obj.predict_state(X_current, acc_meas, gyro_meas, dt);
            
            % 协方差预测 - 基于几何线性化
            P_pred = obj.predict_covariance(X_current, P_current, acc_meas, gyro_meas, dt, Q);
            
            % 确保几何约束
            X_pred = ManifoldOperations.enforce_constraints(X_pred);
            P_pred = obj.ensure_covariance_symmetry(P_pred);
            
            % 记录性能统计
            obj.prediction_times = [obj.prediction_times, toc];
            
            % 约束检查
            if ~ManifoldOperations.is_valid_state(X_pred)
                obj.constraint_violations = obj.constraint_violations + 1;
                warning('预测状态违反约束，已修正');
            end
            
            % 保存调试信息
            if obj.debug_enabled
                obj.last_prediction = struct();
                obj.last_prediction.X_input = X_current;
                obj.last_prediction.X_output = X_pred;
                obj.last_prediction.dt = dt;
                obj.last_prediction.imu = [acc_meas; gyro_meas];
            end
        end
        
        function X_pred = predict_state(obj, X_current, acc_meas, gyro_meas, dt)
            % 状态预测 - 基于几何随机微分方程
            %
            % 实现文章公式(1)-(2)的Munthe-Kaas积分方法
            %
            % 输入:
            %   X_current - 17x1 当前状态向量
            %   acc_meas - 3x1 加速度计测量
            %   gyro_meas - 3x1 陀螺仪测量  
            %   dt - 时间间隔
            %
            % 输出:
            %   X_pred - 17x1 预测状态向量
            
            % 提取当前状态分量
            [pos, vel, att_euler, ba, bg, clk] = ...
                ManifoldOperations.extract_state_components(X_current);
            
            % 当前姿态旋转矩阵
            R_current = SO3.from_euler(att_euler, 'ZYX');
            
            % 根据积分方法选择
            switch lower(obj.integration_method)
                case 'rk4'
                    X_pred = obj.runge_kutta_4th_order(X_current, acc_meas, gyro_meas, dt);
                case 'munthe_kaas'
                    X_pred = obj.munthe_kaas_integration(X_current, acc_meas, gyro_meas, dt);
                case 'euler'
                    X_pred = obj.euler_integration(X_current, acc_meas, gyro_meas, dt);
                otherwise
                    % 默认使用4阶龙格-库塔方法
                    X_pred = obj.runge_kutta_4th_order(X_current, acc_meas, gyro_meas, dt);
            end
        end
        
        function X_next = munthe_kaas_integration(obj, X_current, acc_meas, gyro_meas, dt)
            % Munthe-Kaas积分方法 - 保持李群结构
            %
            %
            % 输入:
            %   X_current - 当前状态
            %   acc_meas, gyro_meas - IMU测量
            %   dt - 时间步长
            %
            % 输出:
            %   X_next - 下一时刻状态
            
            % 计算k1 (初始切向量)
            k1 = obj.compute_drift_vector(X_current, acc_meas, gyro_meas) * dt;
            
            % 计算中点状态
            X_mid = ManifoldOperations.boxplus(X_current, k1/2);
            
            % 计算k2 (中点切向量)  
            k2 = obj.compute_drift_vector(X_mid, acc_meas, gyro_meas) * dt;
            
            % 最终状态更新
            X_next = ManifoldOperations.boxplus(X_current, k2);
        end
        
        function X_next = runge_kutta_4th_order(obj, X_current, acc_meas, gyro_meas, dt)