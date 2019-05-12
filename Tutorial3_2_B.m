clear 
%Part 2a
rate = 20;
dt = 0.1E-3; % binsize
tmax = 100; %max time
tvec = (0:dt:tmax); %time vector

%spikes = alt_poissrnd(tvec);
spikes = rand(1,length(tvec)-1)<rate*dt;

% spike_times = find(spikes)*dt;
% ISI = diff(spike_times);

E_L = -70E-3;
V_th = -50E-3;
V_reset = -80E-3;
Delta_th = 2E-3;
G_L = 10E-9;
C_m = 100E-12;
a = 2E-9;
%b = 0;
b = 1E-9;
tau_SRA = 150E-3;
dt = 0.1E-3;
t_max = 100;

sigma = 50E-12;
std_factor = sigma/sqrt(dt);
t = 0:dt:t_max;
I_app = randn(1,length(t-1))*std_factor;
I_SRA = zeros(size(t));
AP_cnt = 0;
trial = 0;
V_m = zeros(1,length(t)-1);
spike = zeros(1, length(t)-1);

for n = 2:length(t)-1
    trial = trial + 1;
    V_m(n) = ((G_L/C_m) * (E_L - V_m(n-1) + Delta_th*exp((V_m(n-1)-V_th)/Delta_th)) - I_SRA(n-1)/C_m + I_app(n-1)/C_m)*dt + V_m(n-1);
    I_SRA(n) = ((a*(V_m(n-1) - E_L) - I_SRA(n-1))/tau_SRA) * dt + I_SRA(n-1);
    if V_m(n) > V_th
        V_m(n) = V_reset;
        I_SRA(n) = I_SRA(n) + b;
        spike(n) = 1;
        AP_cnt = AP_cnt + 1;
    end 
end


ISI = diff(find(spike)).*dt;
spikenew = expandspikebin(spike, dt,100E-3);
Fano = var(spikenew)./mean(spikenew);
cv = std(ISI)/mean(ISI);

figure(1)
histogram(ISI,25);
xlabel('ISI(sec)')
ylabel('frequency')
