%time vector
dt = 0.1E-3;t_max = 1.5;
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

for n = 2:length(t)
   %if  (n/length(t))> (1/3) && (n/length(t))<(2/3)
    if  (n/length(t))> (1/5) && (n/length(t))<(4/5)

        I_0 = 500E-12;
    else
        I_0 = 0;
    end
    I_app(n) = I_0;
    V_m(n) = ((G_L/C_m) * (E_L - V_m(n-1) + Delta_th*exp((V_m(n-1)-V_th)/Delta_th)) - I_SRA(n-1)/C_m + I_app(n-1)/C_m)*dt + V_m(n-1);
    %V_m(n) = (G_L/C_m)*dt*(E_L - V_m(n-1) + Delta_th*exp((V_m(n-1)-V_th)/Delta_th))-(I_SRA(n-1)*dt)/C_m + (I_app(n-1)*dt)/C_m + V_m(n-1);
    I_SRA(n) = ((a*(V_m(n-1) - E_L) - I_SRA(n-1))/tau_SRA) * dt + I_SRA(n-1);
    if V_m(n) > V_max
        V_m(n) = V_reset;
        I_SRA(n) = I_SRA(n) + b;
    end
end


subplot(2,1,1);
plot(t,I_app);
%xlabel('Time/(s)');
ylabel('Current')

subplot(2,1,2);
plot(t,V_m);
xlabel('Time/(s)');
ylabel('Membrane Potential');






