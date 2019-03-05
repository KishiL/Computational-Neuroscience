clear
%Parameters
E_L = -75E-3;
V_th = -50E-3;
V_reset = -80E-3;
R_m = 100E6;
C_m = 100E-12;
E_K = -80E-3;
DeltaG_sra = 1E-9;
tau_sra = 200E-3;

%time vector
dt = 0.1E-3;t_max = 5;
t = [0:dt:t_max];
DeltaI = 15E-12;

V = zeros(size(t));
V(1) = E_L;
G_sra = zeros(size(t));
G_sra(1) = 0;
FR_ss = zeros(1,20);
I_app = ones(size(t));

%sz = size(I_app(1:2500));
% for n = 1:1:21
%     a = 1+2500*(n-1);
%     b = 2500*n;
%     I_0 = 500E-12 + (n-1)*DeltaI;
%     I_app(a:b) = ones(sz)*I_0;
% end

trial = 0;

for I_0 = 250E-12:DeltaI:550E-12
    trial = trial + 1;
    V = zeros(size(t));
    ap = zeros(size(t));
    for n = 2:length(t)
        %I_0 = 500E-12;
        I_app = ones(size(t))*I_0;
        AP_cnt = 0;
        V(n) = V(n-1) + ((E_L - V(n-1))/(R_m*C_m) + G_sra(n-1)*(E_K-V(n-1))/C_m + I_app(n-1)/C_m)*dt;
        G_sra(n) = (-G_sra(n-1)/tau_sra)*dt + G_sra(n-1);
        if V(n) > V_th
            V(n) = V_reset;
            G_sra(n) = G_sra(n) + DeltaG_sra;
            ap(n) = 1;
            AP_cnt = AP_cnt + 1;
        end
        
    end
    I_in(trial) = I_0;
    FR_ss(trial) = sum(ap(30001:50001))/2;
    ISI_init = diff(find(ap,2)) * dt;
    if ~isnan(ISI_init)
        FR_init(trial) = 1/ISI_init;
    end
end

figure(1)
hold off
plot(I_in,FR_ss);
hold on
stem(I_in,FR_init)
legend('SS rate','initial rate','Location','Northwest');
xlabel('Input Current');
ylabel('Firing Rate')

