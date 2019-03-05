%time vector
dt = 5E-7;t_max = 1;
t = [0:dt:t_max];
DeltaI = 15E-12;

%Parameters
E_L = -75E-3;
E_K = -80E-3;
V_reset = -60E-3;
V_max = 50E-3;
Delta_th = 10E-3;
G_L = 10E-9;
C_m = 100E-12;
% a = 2E-9;
% b = 0.02E-9;
tau_SRA = 200E-3;

G_ref = zeros(size(t));
G_ref(1) = 0;
DeltaG = 60E-9;
tau_Gref = 0.5E-3;

I_app = zeros(size(t));
tau_th = 2E-3;
trial = 0;

Vmax_th = 150E-3;

for I_0 = 150E-12:DeltaI:550E-12
    trial = trial + 1;
    V_m = zeros(size(t));
    V_m(1) = E_L;
    V_th = zeros(size(t));
    V_th(1) = -50E-3;
    ap = zeros(size(t));
    I_app = ones(size(t))*I_0;
    
    for n= 2:length(t)
        AP_cnt = 0;
        V_th(n) = V_th(n-1) + ((V_th(1)-V_th(n-1))/tau_th)*dt;
        G_ref(n) = G_ref(n-1) + (-G_ref(n-1)/tau_Gref)*dt;
        V_m(n) = ((G_L/C_m) * (E_L - V_m(n-1) + Delta_th*exp((V_m(n-1)-V_th(n-1))/Delta_th))+ (G_ref(n-1)*(E_K-V_m(n-1)))/C_m + I_app(n-1)/C_m)*dt + V_m(n-1);
        if V_m(n) > Vmax_th + V_th(n)
            V_m(n) = V_reset;
            G_ref(n) = G_ref(n) + DeltaG;
            V_th(n) = 4*Vmax_th;
            ap(n) = 1;
            AP_cnt = AP_cnt + 1;
        end
    end
    I_in(trial) = I_0;
    ISI = diff(find(ap)) * dt;
    if ~isempty(ISI)
        FR_ss(trial) = 1/ISI(end);
    end
    meanV_m(trial) = mean(V_m);
    %     if ~isnan(ISI_init)
    %         FR_init(trial) = 1/ISI_init;
    %     end
    
    %FR_ss(trial) = 1/FR_ss;
    
    disp(trial)
end


subplot(2,2,1);
plot(t(1:500000),1E3*V_m(1:500000));
xlabel('Time/(s)');
ylabel('Membrane Potential(mV)');

subplot(2,2,2);
plot(I_in,FR_ss);
xlabel('Input Current');
ylabel('Firing Rate(Hz)');

subplot(2,2,3);
plot(FR_ss,1E3*meanV_m);
xlabel('Firing Rate(Hz)');
ylabel('Mean Membrane Potential(mV)');


subplot(2,2,4);
plot(I_in,1E3*meanV_m);
xlabel('Input Current');
ylabel('Mean Membrane Potential(mV)')
