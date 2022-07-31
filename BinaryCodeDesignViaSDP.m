clear
close all

%% 基于最小化峰值非周期自相关旁瓣电平（PSL）的恒模波形优化程序
% s = x + i*y, s为N维复数向量, x为N维实数向量，y为N维实数向量
% N   order    PSL      COPT5.02cputime(s)  MOSEK10.0cputime(s) 
% 5     3    0.7071       225.70               200.27
% 4     3    0.5          20.72                22.28
% 3     3    0            1.20                 2.17

N = 3; % 码元数
x = sdpvar(N,1);  % 实数部分
y = sdpvar(N,1);  % 虚数部分 
h = sdpvar;  %目标函数

F = []; % 约束条件
for n = 1:N
    F = [F, x(n)^2 + y(n)^2 == 1]; % 恒模约束条件  |s(n)| == 1, n=1,2,...,N
end 
for k=N-2:-1:1 %非周期自相关函数旁瓣电平极小化约束
    g = 0;
    gi = 0;
    for i=1:N-k
        g=g+x(i)*x(i+k)+y(i)*y(i+k);
        gi=gi+x(i)*y(i+k)-y(i)*x(i+k);
    end
    F = [F, g^2+gi^2<=h]; % 非周期自相关函数旁瓣电平都小于h
end

order = 3;
opts = sdpsettings('solver','mosek','verbose',2);
% opts = sdpsettings('solver','copt','verbose',2, 'copt.SDPMethod', 2);
opts = sdpsettings('solver','copt','verbose',2);
solvemoment(F, h, opts, order); %极小化h

PSL = (value(h))^0.5

% s = double(x) + 1i*double(y);
% Y = xcorr(s);
% plot(abs(Y))
% hold on
% plot(abs(s))
