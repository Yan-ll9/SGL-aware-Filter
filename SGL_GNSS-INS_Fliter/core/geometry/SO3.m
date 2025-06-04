classdef SO3 < handle
% SO3 特殊正交群 SO(3) 的MATLAB实现
%
% 实现SO(3)上的所有几何操作，包括：
% - 指数映射和对数映射 (Exp, Log)
% - 流形操作符 (boxplus, boxminus, oplus)
% - Rodrigues公式
% - 几何雅可比矩阵
% - 随机元素生成
%
% 
%
% 作者: Y
% 日期: 2025.6

    methods (Static)
        
        function R = exp(omega)
            % EXP SO(3)指数映射: so(3) -> SO(3)
            % 
            % 输入:
            %   omega - 3x1向量，李代数so(3)中的元素（角速度向量）
            %
            % 输出:
            %   R - 3x3旋转矩阵，SO(3)中的元素
            %
            % 实现Rodrigues公式（文章公式3）
            
            if nargin == 0
                R = eye(3);
                return;
            end
            
            omega = omega(:); % 确保是列向量
            
            % 计算角度的模长
            angle = norm(omega);
            
            if angle < 1e-8
                % 小角度近似，避免数值问题
                R = eye(3) + SO3.skew(omega);
                return;
            end
            
            % Rodrigues公式实现（文章公式3）
            axis = omega / angle;
            K = SO3.skew(axis);
            
            R = eye(3) + sin(angle) * K + (1 - cos(angle)) * K^2;
        end
        
        function omega = log(R)
            % LOG SO(3)对数映射: SO(3) -> so(3)
            %
            % 输入:
            %   R - 3x3旋转矩阵
            %
            % 输出:
            %   omega - 3x1向量，对应的李代数元素
            
            % 检查输入有效性
            if ~SO3.is_valid_rotation(R)
                warning('输入不是有效的旋转矩阵');
                R = SO3.nearest_rotation(R);
            end
            
            % 计算旋转角度
            trace_R = trace(R);
            % 防止数值误差导致的trace超出[-1,3]范围
            trace_R = max(-1, min(3, trace_R));
            
            angle = acos((trace_R - 1) / 2);
            
            if angle < 1e-8
                % 小角度情况
                omega = SO3.vee((R - R') / 2);
                return;
            end
            
            if abs(angle - pi) < 1e-8
                % 角度接近π的特殊情况
                % 寻找最大的对角元素
                [~, idx] = max(diag(R));
                
                % 构造轴向量
                axis = zeros(3, 1);
                axis(idx) = sqrt((R(idx, idx) + 1) / 2);
                
                % 确定其他分量的符号
                for i = 1:3
                    if i ~= idx
                        axis(i) = R(idx, i) / (2 * axis(idx));
                    end
                end
                
                omega = angle * axis;
                return;
            end
            
            % 一般情况
            axis = SO3.vee(R - R') / (2 * sin(angle));
            omega = angle * axis;
        end
        
        function R_new = boxplus(R, xi)
            % BOXPLUS 流形加法操作符 ⊞
            %
            % 输入:
            %   R - 3x3旋转矩阵（流形上的点）
            %   xi - 3x1向量（切空间中的扰动）
            %
            % 输出:
            %   R_new - 更新后的旋转矩阵
            %
            % 实现: R ⊞ ξ = R · Exp(ξ)
            
            R_new = R * SO3.exp(xi);
        end
        
        function xi = boxminus(R2, R1)
            % BOXMINUS 流形减法操作符 ⊟
            %
            % 输入:
            %   R2, R1 - 3x3旋转矩阵
            %
            % 输出:
            %   xi - 3x1向量，切空间中的差值
            %
            % 实现: R2 ⊟ R1 = Log(R1^T · R2)
            
            xi = SO3.log(R1' * R2);
        end
        
        function R_new = oplus(R, v)
            % OPLUS 外部输入操作符 ⊕
            %
            % 对于SO(3)，⊕操作与⊞操作相同
            %
            % 输入:
            %   R - 3x3旋转矩阵
            %   v - 3x1向量（外部输入）
            %
            % 输出:
            %   R_new - 更新后的旋转矩阵
            
            R_new = SO3.boxplus(R, v);
        end
        
        function K = skew(v)
            % SKEW 反对称矩阵算子 (·)×
            %
            % 输入:
            %   v - 3x1向量
            %
            % 输出:
            %   K - 3x3反对称矩阵
            
            if length(v) ~= 3
                error('输入向量必须是3维');
            end
            
            v = v(:); % 确保是列向量
            
            K = [0,    -v(3),  v(2);
                 v(3),  0,    -v(1);
                -v(2),  v(1),  0];
        end
        
        function v = vee(K)
            % VEE 反对称矩阵的逆算子 (·)∨
            %
            % 输入:
            %   K - 3x3反对称矩阵
            %
            % 输出:
            %   v - 3x1向量
            
            if size(K, 1) ~= 3 || size(K, 2) ~= 3
                error('输入必须是3x3矩阵');
            end
            
            v = [K(3, 2); K(1, 3); K(2, 1)];
        end
        
        function J = left_jacobian(omega)
            % LEFT_JACOBIAN 左雅可比矩阵
            %
            % 输入:
            %   omega - 3x1向量
            %
            % 输出:
            %   J - 3x3左雅可比矩阵
            %
            % 用于误差状态线性化
            
            omega = omega(:);
            angle = norm(omega);
            
            if angle < 1e-8
                J = eye(3) + 0.5 * SO3.skew(omega);
                return;
            end
            
            axis = omega / angle;
            s = sin(angle);
            c = cos(angle);
            
            J = (s / angle) * eye(3) + ...
                (1 - s / angle) * (axis * axis') + ...
                ((1 - c) / angle) * SO3.skew(axis);
        end
        
        function J_inv = left_jacobian_inv(omega)
            % LEFT_JACOBIAN_INV 左雅可比矩阵的逆
            %
            % 输入:
            %   omega - 3x1向量
            %
            % 输出:
            %   J_inv - 3x3左雅可比矩阵的逆
            
            omega = omega(:);
            angle = norm(omega);
            
            if angle < 1e-8
                J_inv = eye(3) - 0.5 * SO3.skew(omega);
                return;
            end
            
            axis = omega / angle;
            half_angle = 0.5 * angle;
            cot_half = cot(half_angle);
            
            J_inv = half_angle * cot_half * eye(3) + ...
                    (1 - half_angle * cot_half) * (axis * axis') - ...
                    half_angle * SO3.skew(axis);
        end
        
        function Ad = adjoint(R)
            % ADJOINT 伴随表示 Ad_R
            %
            % 输入:
            %   R - 3x3旋转矩阵
            %
            % 输出:
            %   Ad - 3x3伴随矩阵
            %
            % 对于SO(3): Ad_R = R
            
            Ad = R;
        end
        
        function ad = adjoint_algebra(omega)
            % ADJOINT_ALGEBRA 李代数的伴随表示 ad_ω
            %
            % 输入:
            %   omega - 3x1向量
            %
            % 输出:
            %   ad - 3x3伴随矩阵
            %
            % 对于so(3): ad_ω = ω×
            
            ad = SO3.skew(omega);
        end
        
        function R = random()
            % RANDOM 生成随机旋转矩阵
            %
            % 输出:
            %   R - 3x3随机旋转矩阵
            %
            % 使用均匀分布的方法
            
            % 生成随机四元数
            u = rand(3, 1);
            
            q = [sqrt(1 - u(1)) * sin(2 * pi * u(2));
                 sqrt(1 - u(1)) * cos(2 * pi * u(2));
                 sqrt(u(1)) * sin(2 * pi * u(3));
                 sqrt(u(1)) * cos(2 * pi * u(3))];
            
            % 四元数转旋转矩阵
            R = SO3.quat2rotm(q);
        end
        
        function R = from_euler(euler, seq)
            % FROM_EULER 从欧拉角创建旋转矩阵
            %
            % 输入:
            %   euler - 3x1向量，欧拉角 [rad]
            %   seq - 字符串，欧拉角序列（默认'ZYX'）
            %
            % 输出:
            %   R - 3x3旋转矩阵
            
            if nargin < 2
                seq = 'ZYX';
            end
            
            phi = euler(1);   % 绕Z轴
            theta = euler(2); % 绕Y轴  
            psi = euler(3);   % 绕X轴
            
            switch upper(seq)
                case 'ZYX'
                    Rz = [cos(phi), -sin(phi), 0;
                          sin(phi),  cos(phi), 0;
                          0,         0,        1];
                    
                    Ry = [cos(theta),  0, sin(theta);
                          0,           1, 0;
                         -sin(theta),  0, cos(theta)];
                    
                    Rx = [1, 0,         0;
                          0, cos(psi), -sin(psi);
                          0, sin(psi),  cos(psi)];
                    
                    R = Rz * Ry * Rx;
                    
                otherwise
                    error('不支持的欧拉角序列: %s', seq);
            end
        end
        
        function euler = to_euler(R, seq)
            % TO_EULER 旋转矩阵转欧拉角
            %
            % 输入:
            %   R - 3x3旋转矩阵
            %   seq - 字符串，欧拉角序列（默认'ZYX'）
            %
            % 输出:
            %   euler - 3x1向量，欧拉角 [rad]
            
            if nargin < 2
                seq = 'ZYX';
            end
            
            switch upper(seq)
                case 'ZYX'
                    % ZYX欧拉角（Roll-Pitch-Yaw）
                    sy = sqrt(R(1,1)^2 + R(2,1)^2);
                    
                    singular = sy < 1e-6;
                    
                    if ~singular
                        x = atan2(R(3,2), R(3,3));
                        y = atan2(-R(3,1), sy);
                        z = atan2(R(2,1), R(1,1));
                    else
                        x = atan2(-R(2,3), R(2,2));
                        y = atan2(-R(3,1), sy);
                        z = 0;
                    end
                    
                    euler = [z; y; x]; % [Yaw; Pitch; Roll]
                    
                otherwise
                    error('不支持的欧拉角序列: %s', seq);
            end
        end
        
        function R = from_axis_angle(axis, angle)
            % FROM_AXIS_ANGLE 从轴角表示创建旋转矩阵
            %
            % 输入:
            %   axis - 3x1单位向量，旋转轴
            %   angle - 标量，旋转角度 [rad]
            %
            % 输出:
            %   R - 3x3旋转矩阵
            
            axis = axis(:) / norm(axis); % 归一化
            omega = angle * axis;
            R = SO3.exp(omega);
        end
        
        function [axis, angle] = to_axis_angle(R)
            % TO_AXIS_ANGLE 旋转矩阵转轴角表示
            %
            % 输入:
            %   R - 3x3旋转矩阵
            %
            % 输出:
            %   axis - 3x1单位向量，旋转轴
            %   angle - 标量，旋转角度 [rad]
            
            omega = SO3.log(R);
            angle = norm(omega);
            
            if angle < 1e-8
                axis = [1; 0; 0]; % 任意轴
                angle = 0;
            else
                axis = omega / angle;
            end
        end
        
        function R = quat2rotm(q)
            % QUAT2ROTM 四元数转旋转矩阵
            %
            % 输入:
            %   q - 4x1四元数 [qx; qy; qz; qw]
            %
            % 输出:
            %   R - 3x3旋转矩阵
            
            q = q(:) / norm(q); % 归一化
            
            qx = q(1); qy = q(2); qz = q(3); qw = q(4);
            
            R = [1 - 2*(qy^2 + qz^2), 2*(qx*qy - qz*qw), 2*(qx*qz + qy*qw);
                 2*(qx*qy + qz*qw), 1 - 2*(qx^2 + qz^2), 2*(qy*qz - qx*qw);
                 2*(qx*qz - qy*qw), 2*(qy*qz + qx*qw), 1 - 2*(qx^2 + qy^2)];
        end
        
        function q = rotm2quat(R)
            % ROTM2QUAT 旋转矩阵转四元数
            %
            % 输入:
            %   R - 3x3旋转矩阵
            %
            % 输出:
            %   q - 4x1四元数 [qx; qy; qz; qw]
            
            trace_R = trace(R);
            
            if trace_R > 0
                s = sqrt(trace_R + 1) * 2; % s = 4 * qw
                qw = 0.25 * s;
                qx = (R(3,2) - R(2,3)) / s;
                qy = (R(1,3) - R(3,1)) / s;
                qz = (R(2,1) - R(1,2)) / s;
            elseif R(1,1) > R(2,2) && R(1,1) > R(3,3)
                s = sqrt(1 + R(1,1) - R(2,2) - R(3,3)) * 2; % s = 4 * qx
                qw = (R(3,2) - R(2,3)) / s;
                qx = 0.25 * s;
                qy = (R(1,2) + R(2,1)) / s;
                qz = (R(1,3) + R(3,1)) / s;
            elseif R(2,2) > R(3,3)
                s = sqrt(1 + R(2,2) - R(1,1) - R(3,3)) * 2; % s = 4 * qy
                qw = (R(1,3) - R(3,1)) / s;
                qx = (R(1,2) + R(2,1)) / s;
                qy = 0.25 * s;
                qz = (R(2,3) + R(3,2)) / s;
            else
                s = sqrt(1 + R(3,3) - R(1,1) - R(2,2)) * 2; % s = 4 * qz
                qw = (R(2,1) - R(1,2)) / s;
                qx = (R(1,3) + R(3,1)) / s;
                qy = (R(2,3) + R(3,2)) / s;
                qz = 0.25 * s;
            end
            
            q = [qx; qy; qz; qw];
        end
        
        function valid = is_valid_rotation(R)
            % IS_VALID_ROTATION 检查是否为有效旋转矩阵
            %
            % 输入:
            %   R - 3x3矩阵
            %
            % 输出:
            %   valid - 逻辑值，是否为有效旋转矩阵
            
            if size(R, 1) ~= 3 || size(R, 2) ~= 3
                valid = false;
                return;
            end
            
            % 检查正交性: R'*R = I
            orth_error = norm(R' * R - eye(3), 'fro');
            
            % 检查行列式: det(R) = 1
            det_error = abs(det(R) - 1);
            
            valid = (orth_error < 1e-6) && (det_error < 1e-6);
        end
        
        function R_nearest = nearest_rotation(M)
            % NEAREST_ROTATION 寻找最近的旋转矩阵
            %
            % 输入:
            %   M - 3x3矩阵
            %
            % 输出:
            %   R_nearest - 最近的旋转矩阵（Frobenius范数意义下）
            %
            % 使用SVD分解方法
            
            [U, ~, V] = svd(M);
            R_nearest = U * V';
            
            % 确保行列式为1（右手坐标系）
            if det(R_nearest) < 0
                V(:, 3) = -V(:, 3);
                R_nearest = U * V';
            end
        end
        
        function dist = distance(R1, R2)
            % DISTANCE 计算两个旋转矩阵之间的距离
            %
            % 输入:
            %   R1, R2 - 3x3旋转矩阵
            %
            % 输出:
            %   dist - 测地距离 [rad]
            
            omega_diff = SO3.log(R1' * R2);
            dist = norm(omega_diff);
        end
        
        function omega_dot = kinematics(R, omega_body)
            % KINEMATICS SO(3)运动学方程
            %
            % 输入:
            %   R - 当前旋转矩阵
            %   omega_body - 体坐标系角速度
            %
            % 输出:
            %   omega_dot - 李代数空间的导数
            %
            % 用于积分更新: ω̇ = J_l^(-1) * ω_body
            
            % 当前状态的李代数表示
            omega_current = SO3.log(R);
            
            % 左雅可比矩阵的逆
            J_l_inv = SO3.left_jacobian_inv(omega_current);
            
            % 运动学方程
            omega_dot = J_l_inv * omega_body;
        end
        
        function display_matrix(R, name)
            % DISPLAY_MATRIX 显示旋转矩阵（调试用）
            %
            % 输入:
            %   R - 3x3旋转矩阵
            %   name - 矩阵名称（可选）
            
            if nargin < 2
                name = 'Rotation Matrix';
            end
            
            fprintf('\n%s:\n', name);
            fprintf('[%8.4f %8.4f %8.4f]\n', R(1,:));
            fprintf('[%8.4f %8.4f %8.4f]\n', R(2,:));
            fprintf('[%8.4f %8.4f %8.4f]\n', R(3,:));
            
            % 验证正交性和行列式
            orth_error = norm(R' * R - eye(3), 'fro');
            det_val = det(R);
            
            fprintf('正交性误差: %.2e\n', orth_error);
            fprintf('行列式: %.6f\n', det_val);
            fprintf('有效性: %s\n', logical_to_string(SO3.is_valid_rotation(R)));
        end
        
    end
    
    methods (Static, Hidden)
        % 隐藏的辅助方法
        
        function run_tests()
            % RUN_TESTS 运行SO(3)类的单元测试
            
            fprintf('=== SO(3)类单元测试 ===\n');
            
            % 测试1: 指数/对数映射的逆性
            fprintf('测试1: 指数/对数映射逆性...');
            omega_test = [0.1; 0.2; 0.3];
            R_test = SO3.exp(omega_test);
            omega_recovered = SO3.log(R_test);
            error1 = norm(omega_test - omega_recovered);
            
            if error1 < 1e-10
                fprintf(' 通过 (误差: %.2e)\n', error1);
            else
                fprintf(' 失败 (误差: %.2e)\n', error1);
            end
            
            % 测试2: 流形操作符一致性
            fprintf('测试2: 流形操作符一致性...');
            R1 = SO3.random();
            R2 = SO3.random();
            xi = SO3.boxminus(R2, R1);
            R2_recovered = SO3.boxplus(R1, xi);
            error2 = norm(R2 - R2_recovered, 'fro');
            
            if error2 < 1e-10
                fprintf(' 通过 (误差: %.2e)\n', error2);
            else
                fprintf(' 失败 (误差: %.2e)\n', error2);
            end
            
            % 测试3: 雅可比矩阵性质
            fprintf('测试3: 左雅可比矩阵逆性...');
            omega_test = [0.05; 0.1; 0.15];
            J_l = SO3.left_jacobian(omega_test);
            J_l_inv = SO3.left_jacobian_inv(omega_test);
            error3 = norm(J_l * J_l_inv - eye(3), 'fro');
            
            if error3 < 1e-10
                fprintf(' 通过 (误差: %.2e)\n', error3);
            else
                fprintf(' 失败 (误差: %.2e)\n', error3);
            end
            
            % 测试4: 四元数转换
            fprintf('测试4: 四元数转换一致性...');
            R_test = SO3.random();
            q = SO3.rotm2quat(R_test);
            R_recovered = SO3.quat2rotm(q);
            error4 = norm(R_test - R_recovered, 'fro');
            
            if error4 < 1e-10
                fprintf(' 通过 (误差: %.2e)\n', error4);
            else
                fprintf(' 失败 (误差: %.2e)\n', error4);
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