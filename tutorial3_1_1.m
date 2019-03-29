E_L = -60E-3;
V_th = -50E-3;
V_reset = -80E-3;
V_max = 50E-3;
Delta_th = 2E-3;
G_L = 8E-9;
C_m = 100E-12;
a = 10E-9;
b = 0.5E-9;
tau_SRA = 50E-3;

%time vector
dt = 0.02E-3;
t_max = 40000*5E-3;
t = [0:dt:t_max];
tsteps = 250;

I_rand = (rand(1,40000)-0.5)*1E-8;
I_app = ones(size(t));

for k = 1:length(I_rand)
    I_app(((k-1)*tsteps+1):(k*tsteps)) = ones(1,tsteps)*I_rand(k);
end

% plot(t,I_app);
% axis([0 100 -5E-10 5E-10]);
I_SRA = zeros(size(t));

AP_cnt = 0;
trial = 0;
V_m = zeros(1,length(t)-1);
spike = zeros(1, length(t)-1);

for n = 2:length(t)-1
    trial = trial + 1;
    V_m(n) = ((G_L/C_m) * (E_L - V_m(n-1) + Delta_th*exp((V_m(n-1)-V_th)/Delta_th)) - I_SRA(n-1)/C_m + I_app(n-1)/C_m)*dt + V_m(n-1);
    I_SRA(n) = ((a*(V_m(n-1) - E_L) - I_SRA(n-1))/tau_SRA) * dt + I_SRA(n-1);
    if V_m(n) > V_max
        V_m(n) = V_reset;
        I_SRA(n) = I_SRA(n) + b;
        spike(n) = 1;
        AP_cnt = AP_cnt + 1;
    end 
end

newI_app = expandbin(I_app,0.02E-3,1E-3);
newspike = expandbin(spike,0.02E-3,1E-3);
[sta_v,tcorr_v] = STA(newI_app,newspike,1E-3);
figure(1);
plot(-tcorr_v,sta_v);
xlabel('tcorr');
ylabel('STA');


