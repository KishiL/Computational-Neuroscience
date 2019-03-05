%time vector
dt = 0.1E-3;t_max = 5;
t = [0:dt:t_max];
DeltaI = 15E-12;

%Parameters
E_L = -75E-3;
V_th = -50E-3;
V_reset = -80E-3;
V_max = 50E-3;
Delta_th = 2E-3;
G_L = 10E-9;
C_m = 100E-12;
a = 2E-9;
b = 0.02E-9;
tau_SRA = 200E-3;

V_m = zeros(size(t));
V_m(1) = E_L;
I_SRA = zeros(size(t));
I_SRA(1) = 0;
I_app = zeros(size(t));
trial = 0;

for I_0 = 250E-12:DeltaI:550E-12
    trial = trial + 1;
    V_m = zeros(size(t));
    ap = zeros(size(t));
    for n = 2:length(t)
        I_app = ones(size(t))*I_0;
        AP_cnt = 0;
        V_m(n) = ((G_L/C_m) * (E_L - V_m(n-1) + Delta_th*exp((V_m(n-1)-V_th)/Delta_th)) - I_SRA(n-1)/C_m + I_app(n-1)/C_m)*dt + V_m(n-1);
        I_SRA(n) = ((a*(V_m(n-1) - E_L) - I_SRA(n-1))/tau_SRA) * dt + I_SRA(n-1);
        if V_m(n) > V_max
            V_m(n) = V_reset;
            I_SRA(n) = I_SRA(n) + b;
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
xlabel('Input current');
ylabel('Firing rate')



