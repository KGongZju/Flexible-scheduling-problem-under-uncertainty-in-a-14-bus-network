%% 14节点网络，测试不确定性影响下的灵活性调度问题
%% 案例设计
% 标准测试网络参数
clear;clc;
mpc = case14;
Nbus = size(mpc.bus(:,1),1);
Nbranch = size(mpc.branch(:,1),1);
Ngen = size(mpc.gen(:,1),1);
LineI = mpc.branch(:,1);
LingJ = mpc.branch(:,2);
% 分布式资源配置位置
L_es = [0;0;1;0;1;0;0;0;1;0;0;0;1;1] ;
L_pv = [0;1;1;0;0;1;1;1;0;0;0;0;0;0] ;
L_wd = [0;1;0;1;1;0;0;1;1;0;0;0;0;0];
L_dg = [1;1;1;0;0;1;0;1;0;0;0;0;0;0];
P_dg = [332.4;140;100;0;0;100;0;100;0;0;0;0;0;0];
ramp_dg = 20;
% 风光数据
load winddata.mat
load pvData.mat
wind = Winddata(1:24)/100 ;
pv = pvData(1:24)/10;
% 负荷变化数据
load_p = [0.8;0.82;0.83;0.9;0.93;0.95;1;1.02;1.04;1.06;1.1;1.12;1.15;1.12;1.08;1.05;1.03;1.0;1.02;1.04;1.06;1.02;0.94;0.88];
load_ex = 1; %负荷扩大倍数
ld = (mpc.bus(:,3) * load_p')';
% 储能设备
es_pmax = 20;
es_cmax = 100;
es_effi = 0.9;
es_cini = 40;
% 成本系数
c_es = 13.2;
c_dg = [43.03;250;10;0;0;10;0;10;0;0;0;0;0;0]/10;
c_pv1 = 10.1;
c_pv2 = 1.4;
c_wpen = 500;
c_res = 30;
% 时间
time =24;
% 不确定性因素
epsilon = [0.05;0.05;0.05;0.05;0.05;0.05;0.05;0.05;0.05;0.05;0.05;0.05;0.05;0.05];
delta = [0.05;0.05;0.05;0.05;0.05;0.05;0.05;0.05;0.05;0.05;0.05;0.05;0.05;0.05]+0.16;
lmd = [1;1;1;1;1;1;1;1;1;1;1;1;1;1];
%% 决策变量声明
es_ch = sdpvar(time,Nbus); %储能充电
es_dh = sdpvar(time,Nbus); %储能放电
es_cp = sdpvar(time,Nbus); %储能容量
an_u1 = binvar(time,Nbus); %辅助变量-储能系统
dc_dg = sdpvar(time,Nbus); %分布式燃机
dc_pv = sdpvar(time,Nbus); %分布式可控光伏
re_up = sdpvar(time,1); %外购向上容量
re_dn = sdpvar(time,1); %外购向下容量
pen_w = sdpvar(time,Nbus); %惩罚弃风量
dc_wd = sdpvar(time,Nbus); %并网风电量
an_u2 = binvar(time,Nbus); %辅助变量-风电
dc_br = sdpvar(time,Nbranch); %线路潮流
theta = sdpvar(time,Nbus); % 节点相角
ramp_dg_up = sdpvar(time,Nbus); %燃机爬坡
ramp_dg_dn = sdpvar(time,Nbus); %燃机滑坡
%% 模型构建
con = [];
for i = 1:time
    for j = 1:Nbranch
        con = [con, dc_br(i,j) == (theta(i,mpc.branch(j,1))-theta(i,mpc.branch(j,2)))/mpc.branch(j,4)];
    end
end

for i = 1:time
    for j = 1:Nbus
        corr = SearchNodeConnection(LineI, LingJ, j); %与j节点关联的节点，位于第二列
        Output=FindBranchNumber(corr,mpc.branch(:,1:2));
        con = [con, ld(i,j)*load_ex-dc_wd(i,j) <= dc_pv(i,j) + dc_dg(i,j) + es_dh(i,j) - es_ch(i,j) - sum(Output(:,1)'.*dc_br(i,Output(:,2)))];
        con = [con, dc_pv(i,j) + dc_dg(i,j) + es_dh(i,j) - es_ch(i,j) - sum(Output(:,1)'.*dc_br(i,Output(:,2))) <= (ld(i,j)*load_ex-dc_wd(i,j)) * (1 - epsilon(j) * lmd(j)) + delta(j) * abs(ld(i,j)*load_ex-dc_wd(i,j))];
    end
end

for i = 1:time
    con = [con, sum(ld(i,:)*load_ex-dc_wd(i,:)) <= sum(dc_pv(i,:)) + sum(dc_dg(i,:)) + sum(es_dh(i,:)) - sum(es_ch(i,:)) - sum(dc_br(i,:))];
    con = [con, sum(dc_pv(i,:)) + sum(dc_dg(i,:)) + sum(es_dh(i,:)) - sum(es_ch(i,:)) - sum(dc_br(i,:)) <= sum(ld(i,:)*load_ex-dc_wd(i,:) .* (1 - epsilon(j) * lmd(:))') + sum(abs(ld(i,j)*load_ex-dc_wd(i,:)).*delta(:)')];
end

for i = 1:time
    for j = 1:Nbus
        con = [con, 0 <= dc_dg(i,j) + ramp_dg_up(i,j) - ramp_dg_dn(i,j) <= L_dg(j)*P_dg(j)];
        con = [con, 0 <= dc_dg(i,j) <= L_dg(j)*P_dg(j)];
        con = [con, 0<= ramp_dg_up(i,j) <= L_dg(j) * ramp_dg];
        con = [con, 0<= ramp_dg_dn(i,j) <= L_dg(j) * ramp_dg];
        con = [con, 0 <= dc_pv(i,j) <= L_pv(j)*pv(i)];
        con = [con, dc_wd(i,j) == an_u2(i,j) * L_wd(j) * wind(i)];
        con = [con, pen_w(i,j) == (1-an_u2(i,j)) * L_wd(j) * wind(i)];
        con = [con, 0<= es_dh(i,j) <= an_u1(i,j) * L_es(j) * es_pmax];
        con = [con, 0<= es_ch(i,j) <= (1-an_u1(i,j)) * L_es(j) * es_pmax];
        if i == 1
            con = [con, (-1+epsilon(j)*lmd(j))*es_effi * es_ch(i,j) + (1+epsilon(j)*lmd(j)) * es_dh(i,j)/es_effi <= (1) * es_cini];
            con = [con, (1+epsilon(j)*lmd(j))*es_effi * es_ch(i,j) + (-1+epsilon(j)*lmd(j)) * es_dh(i,j)/es_effi <= (1) * (es_cmax - es_cini)];
        else
            con = [con, (-1+epsilon(j)*lmd(j))*es_effi * es_ch(i,j) + (1+epsilon(j)*lmd(j)) * es_dh(i,j)/es_effi <= (1) * es_cp(i-1,j)];
            con = [con, (1+epsilon(j)*lmd(j))*es_effi * es_ch(i,j) + (-1+epsilon(j)*lmd(j)) * es_dh(i,j)/es_effi <= (1) * (es_cmax - es_cp(i-1,j))];
        end
    end
end

for i = 1:time
    con = [con, sum(ramp_dg_up(i,:)) + sum(L_es(:)'*es_pmax + es_ch(i,:) - es_dh(i,:)) + sum(L_pv(:)'*pv(i) - dc_pv(i,:)) +  re_up(i) >= 0.1 * sum(ld(i,:)) + 0.3 * sum(L_wd(:)' .* dc_wd(i,:))];
    con = [con, sum(ramp_dg_dn(i,:)) + sum(L_es(:)'*es_pmax - es_ch(i,:) + es_dh(i,:)) + sum(dc_pv(i,:)) +  re_dn(i) >= 0.1 * sum(ld(i,:)) + 0.3 * sum(L_wd(:)' .* dc_wd(i,:))];
    con = [con, re_up(i) >= 0, re_dn(i) >= 0];
end
es_cost = sdpvar(time,Nbus);
dg_cost = sdpvar(time,Nbus);
re_cost = sdpvar(time,1);
% pv_cost = sdpvar(time,Nbus);
pen_cost = sdpvar(time,Nbus);
for i = 1:time
    for j = 1:Nbus
        con = [con, es_cost(i,j) == c_es * (es_dh(i,j) + es_ch(i,j))];
        con = [con, dg_cost(i,j) == c_dg(j) * dc_dg(i,j)];
        con = [con, re_cost(i) == c_res * (re_up(i)+re_dn(i))];
        pv_cost(i,j) = c_pv1 * dc_pv(i,j) + c_pv2 * (L_pv(j)*pv(i)-dc_pv(i,j))^2;
        con= [con, pen_cost(i,j) == c_wpen * pen_w(i,j)];
    end
end
obj = sum(sum(es_cost+dg_cost+pv_cost+pen_cost)) + sum(re_cost);
ops=sdpsettings('solver','gurobi');
% ops.gurobi.InfUnbdInfo=1;
optimize(con,obj,ops);

dc_br = double(dc_br);
dc_dg = double(dc_dg);
dc_pv = double(dc_pv);
dc_wd = double(dc_wd);
dg_cost = double(dg_cost);
dg_cost = sum(sum(dg_cost));
es_ch = double(es_ch);
es_cost = double(es_cost);
es_cost = sum(sum(es_cost));
es_cp = double(es_cp);
es_dh = double(es_dh);
pen_cost = double(pen_cost);
pen_cost = sum(sum(pen_cost));
pen_w = double(pen_w);
pv_cost = double(pv_cost);
pv_cost = sum(sum(pv_cost));
re_cost = double(re_cost);
re_cost = sum(re_cost);
re_dn = double(re_dn);
re_up = double(re_up);
ramp_dg_up = double(ramp_dg_up);
ramp_dg_dn = double(ramp_dg_dn);
obj = double(obj);

curve_wd_inject = sum(dc_wd');
curve_wd_abondon = sum(pen_w');
curve_pv = sum(dc_pv');
curve_dg = sum(dc_dg');
curve_es = sum(es_dh'-es_ch');
curve_ld = sum(ld');
curve_input = curve_es+curve_dg+curve_pv+curve_wd_inject;

figure
plot(curve_ld);
hold on;
plot(sum(L_pv)*pv)
plot(sum(L_wd)*wind)

figure
stairs(curve_dg)
hold on;
stairs(curve_dg + curve_wd_inject)
% stairs(curve_dg + curve_wd_inject + curve_pv)
stairs(curve_dg + curve_wd_inject + curve_pv + curve_es)
stairs(curve_dg + curve_wd_inject + curve_wd_abondon);

curve_windabandonvalue = [230.513;153.684;153.684;153.684;76.723;76.723;76.723;0;0;0;0;0;0;0;0;0;];
curve_windabandonrate = [0.0345;0.023;0.023;0.023;0.0115;0.0115;0.0115;0;0;0;0;0;0;0;0;0;];
curve_reseveloss = [11.3823;32.603;24.0471;16.8508;43.8217;35.2345;26.6473;49.1211;40.6034;32.0856;23.5678;15.0501;9.4813;5.2421;1.003;0;];
curve_losstime = [1;2;2;2;2;2;2;2;2;2;2;2;1;1;1;0;];
deltavalue=[0.06:0.01:0.21]';
figure
subplot(4,1,1)
plot(deltavalue,curve_windabandonvalue);
subplot(4,1,2)
plot(deltavalue,curve_windabandonrate);
subplot(4,1,3)
plot(deltavalue,curve_reseveloss);
subplot(4,1,4)
plot(deltavalue,curve_losstime);

sum(curve_wd_abondon)
sum(curve_wd_abondon)/sum(curve_wd_abondon+curve_wd_inject)
sum(re_dn+re_up)