clear
%
E_L = -70E-3;
V_th = -50E-3;
V_reset = -80E-3;
Delta_th = 2E-3;
G_L = 10E-9;
C_m = 100E-12;
a = 2E-9;
b = 0E-9;
tau_SRA = 150E-3;

sigma = 20E-12;       
sigma_sp = 5E-12;
sigma_sa = 50E-12;
dt = 0.01E-3;

t_max = 0.5;
% t_max = 0.2;
tvec=0:dt:t_max;
std_factor = sigma/sqrt(dt);
std_sp = sigma_sp/sqrt(dt);
std_sa = sigma_sa/sqrt(dt);

for i = 1:1000
%     I_app_sp(i,:) = randn(1,length(tvec)-1)*std_factor + 0.1E-9;
    I_app_sp(i,:) = randn(1,length(tvec)-1)*std_sp + 0.5E-9;
end


for j = 1:1000
%     I_app_sa(j,:) = randn(1,length(tvec)-1)*std_factor;
    I_app_sa(j,:) = randn(1,length(tvec)-1)*std_sa;
end

% I_app_sp = randn(1,length(tvec)-1)*std_factor + 0.1E-9;
% I_app_sa = randn(1,length(tvec)-1)*std_factor;

for k = 1:1000
    I_SRA = zeros(size(tvec));
    AP_cnt = 0;
    trial = 0;
    V_m = zeros(1,length(tvec)-1);
    spike = zeros(1, length(tvec)-1);
    
    for n = 2:length(tvec)-1
        trial = trial + 1;
        V_m(n) = ((G_L/C_m) * (E_L - V_m(n-1) + Delta_th*exp((V_m(n-1)-V_th)/Delta_th)) - I_SRA(n-1)/C_m + I_app_sp(k,n-1)/C_m)*dt + V_m(n-1);
        I_SRA(n) = ((a*(V_m(n-1) - E_L) - I_SRA(n-1))/tau_SRA) * dt + I_SRA(n-1);
        if V_m(n) > V_th
            V_m(n) = V_reset;
            I_SRA(n) = I_SRA(n) + b;
            spike(n) = 1;
            AP_cnt = AP_cnt + 1;
        end
    end
    spike_sp(k,:)=spike; 
end

for s = 1:1000
    I_SRA = zeros(size(tvec));
    AP_cnt = 0;
    trial = 0;
    V_m = zeros(1,length(tvec)-1);
    spike = zeros(1, length(tvec)-1);
    
    for n = 2:length(tvec)-1
        trial = trial + 1;
        V_m(n) = ((G_L/C_m) * (E_L - V_m(n-1) + Delta_th*exp((V_m(n-1)-V_th)/Delta_th)) - I_SRA(n-1)/C_m + I_app_sa(s,n-1)/C_m)*dt + V_m(n-1);
        I_SRA(n) = ((a*(V_m(n-1) - E_L) - I_SRA(n-1))/tau_SRA) * dt + I_SRA(n-1);
        if V_m(n) > V_th
            V_m(n) = V_reset;
            I_SRA(n) = I_SRA(n) + b;
            spike(n) = 1;
            AP_cnt = AP_cnt + 1;
        end
    end
    spike_sa(s,:)=spike; 
end

% % spike_sum_sp = sum(spike_sp,2);
% % spike_mean_sp = mean(spike_sum);
% 
% fr_sp = mean(sum(spike_sp,2))/t_max;
% 
% 
% spike_hist_sp = histcounts(sum(spike_sp,2));
% figure(1);
% stairs(spike_hist_sp);
% 
% frc_sp= 1-(cumsum(spike_hist_sp)./1000);
% % figure(2);
% % stairs(frc_sp);
% 
% fr_sa = mean(sum(spike_sa,2))/t_max;
% spike_hist_sa = histcounts(sum(spike_sa,2));
% % figure(3);
% % stairs(spike_hist_sa);
% 
% frc_sa= 1-(cumsum(spike_hist_sa)./1000);
% % figure(4);
% % stairs(frc_sa);

spike_sum_sp = sum(spike_sp,2);
spike_sum_sa = sum(spike_sa,2);


BMAX = max(max([spike_sum_sp spike_sum_sa]));
[N1,edges1] = histcounts(spike_sum_sp,'BinLimits',[0,BMAX]);
BW = edges1(2)-edges1(1);
[N2,edges2] = histcounts(spike_sum_sa,'BinWidth', BW,'BinLimits',[0,BMAX]);

P1 = cumsum(N1/1000, 'reverse');
P2 = cumsum(N2/1000, 'reverse');

figure(1);
stairs(edges1(1:end-1), P1)
hold on
stairs(edges1(1:end-1), P2)
xlabel('Spike count,\chi')
ylabel('Probability greater')
legend('\alpha(\chi)','\beta(\chi)')

figure(2);
stairs(N1)
hold on
stairs(N2)
xlabel('Spike count, \chi');
ylabel('Number of trials');
legend('\alpha(\chi)','\beta(\chi)')

figure(3)
plot(P2,P1)
xlabel('Probability of false +ve,\beta(\chi)')
ylabel('Probability of false +ve,\alpha(\chi)')

% 
% figure(1);
% stairs(spike_hist_sp);
% hold on;
% stairs(spike_hist_sa);
% xlabel('Spike count, \chi');
% ylabel('Number of trials');
% 
% figure(3);
% stairs(frc_sp);
% hold on;
% stairs(frc_sa);
% xlabel('Spike count,\chi')
% ylabel('Probability greater')
% legend('\alpha(\chi)','\beta(\chi)')



% figure(1);
% stairs(spike_hist_sp);
% figure(2);
% stairs(frc_sp);

% BMAX = max([spike1_cnt spike2_cnt]);
% [N1,edges1] = histcounts(spike1_cnt,'BinLimits',[0,BMAX]);
% BW = edges1(2)-edges1(1);
% [N2,edges2] = histcounts(spike2_cnt,'BinWidth', BW,'BinLimits',[0,BMAX]);
% 
% P1 = cumsum(N1/1000, 'reverse');
% P2 = cumsum(N2/1000, 'reverse');
% stairs(edges1(1:end-1), P1)
% hold on
% stairs(edges1(1:end-1), P2)


