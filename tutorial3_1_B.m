E_L = -60E-3;
V_th = -50E-3;
V_reset = -80E-3;
V_max = 50E-3;
Delta_th = 2E-3;
G_L = 8E-9;
C_m = 100E-12;
a = 40E-9;
b = 1E-9;
tau_SRA = 50E-3;

%time vector
dt = 0.02E-3;
t_max = 40000*5E-3;
t = 0:dt:t_max;
tsteps = 250;
x0 = 20.5;

S = (rand(40,40000)-0.5)*1E-8;
[S_m,S_n] = size(S);
x_max = S_m;
W = zeros(1,x_max);

for x = 1:x_max
    W(x) = cos(4*pi*((x-x0)/x_max))*exp(-16*((x-x0)/x_max).^2);    
end

I_0 = W*S;

for k = 1:length(I_0)
    I_app(((k-1)*tsteps+1):(k*tsteps)) = ones(1,tsteps)*I_0(k);
end

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

I_sp = zeros(S_m,S_n);
I_sp_ds = zeros(S_m,200000);

%upsample S matrix
for n = 1:40
    for k =1:length(I_0)
        I_sp(n,(k-1)*tsteps+1:k*tsteps) = ones(1,tsteps)*S(n,k);
    end
    I_sp_ds(n,:) = expandbin(I_sp(n,:),dt,1E-3);
end

[sta_v,tcorr_v] = STA_spatial(I_sp_ds,newspike,1E-3);
% imagesc(fliplr(sta));

figure(1);
imagesc((sta_v))

