clear;          % 清除工作区变量
clc;            % 清除命令窗口
close all;      % 关闭所有图形窗口

%% 基本电磁参数定义
Frequency = 10e9;       % 工作频率（10 GHz）
Lightspeed = physconst('LightSpeed');  % 光速（物理常数）
Wavelength = Lightspeed/Frequency;    % 计算波长（λ = c/f）
Wavenumber = 2*pi/Wavelength;         % 波数（k = 2π/λ）

%% 阵列参数配置
N = 25;                  % 阵列元素数量
X = (1:N)*Wavelength*0.5;% 元素位置（半波长间距）
I0 =  ones(1,N);         % 均匀幅度分布（所有元素幅度相同）
alpha = zeros(1,N);      % 相位分布（所有元素同相位）

%% 阵列因子采样参数
Ns = 500;                % 角度采样点数
theta = linspace(-90,90,Ns);  % 角度范围（-90°到90°）

%% 均匀线阵计算
E0 = zeros(1,Ns);        % 预分配存储均匀阵列响应
for num = 1:Ns           % 遍历所有角度点
    % 计算阵列因子：Σ I_n * exp(jkx sinθ + jα_n)
    E0(num) = sum(I0 .* exp(1j*(Wavenumber*X*sind(theta(num)) + alpha))) + 1e-3;
end
E0_dB = db(E0) - max(db(E0));  % 转换为分贝并归一化

%% 切比雪夫线阵计算
E1 = zeros(1,Ns);        % 预分配存储切比雪夫阵列响应
I1 = chebwin(N,30)';     % 生成切比雪夫幅度分布（30 dB旁瓣电平）
for num = 1:Ns
    E1(num) = sum(I1 .* exp(1j*(Wavenumber*X*sind(theta(num)) + alpha))) + 1e-3;
end
E1_dB = db(E1) - max(db(E1));  % 转换为分贝并归一化

%% 泰勒线阵计算
E2 = zeros(1,Ns);        % 预分配存储泰勒阵列响应
I2 = taylorwin(N,4,-30)';  % 生成泰勒幅度分布（4个等旁瓣，-30 dB）
for num = 1:Ns
    E2(num) = sum(I2 .* exp(1j*(Wavenumber*X*sind(theta(num)) + alpha))) + 1e-3;
end
E2_dB = db(E2) - max(db(E2));  % 转换为分贝并归一化

%% 绘制波束图
figure;
plot(theta,E0_dB,'LineWidth',2);  % 可选：绘制均匀阵列
hold on
plot(theta,E1_dB,'LineWidth',2);   % 绘制切雪比夫阵列
hold on
ylim([-40,0]);                    % 设置纵轴范围
grid on
xlabel('\theta(\circ)');          % X轴标签
ylabel('dB');                     % Y轴标签
axis tight;                       % 自动调整坐标轴
figure;
plot(theta,E0_dB,'LineWidth',2);  % 可选：绘制均匀阵列
hold on
plot(theta,E2_dB,'LineWidth',2);   % 绘制泰勒阵列
ylim([-40,0]);                    % 设置纵轴范围
grid on                           % 添加网格
xlabel('\theta(\circ)');          % X轴标签
ylabel('dB');                     % Y轴标签
axis tight;                       % 自动调整坐标轴
figure;
plot(theta,E1_dB,'LineWidth',2);  % 绘制切雪比夫阵列
hold on
plot(theta,E2_dB,'LineWidth',2);   % 绘制泰勒阵列
hold on
ylim([-40,0]);                    % 设置纵轴范围
grid on
legend("切雪比夫阵列","泰勒阵列")
xlabel('\theta(\circ)');          % X轴标签
ylabel('dB');                     % Y轴标签
axis tight;   
%% 绘制幅度分布
figure()
scatter(1:25,I1,'filled');        % 绘制切比雪夫幅度分布
hold on
scatter(1:25,I2,'filled');        % 绘制泰勒幅度分布
legend('Chebyshev (30 dB)','Taylor (4,-30)','Location','best');
xlabel('Element Number');         % X轴标签
ylabel('Amplitude');              % Y轴标签
grid on                           % 添加网格