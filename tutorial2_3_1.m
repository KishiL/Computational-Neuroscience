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
I_0 = 0;

%time vector
dt = 0.1E-3;t_max = 1.5;
t = [0:dt:t_max];

V = zeros(size(t));
V(1) = E_L;
G_sra = zeros(size(t));
G_sra(1) = 0;

%I_0 = 500E-9;
I_app = zeros(size(t));

for n = 2:length(t)
    if  (n/length(t))> (1/3) && (n/length(t))<(2/3)
        I_0 = 500E-12;
    else
        I_0 = 0;
    end
    I_app(n) = I_0;
    V(n) = V(n-1) + ((E_L - V(n-1))/(R_m*C_m) + G_sra(n-1)*(E_K-V(n-1))/C_m + I_app(n-1)/C_m)*dt;
    G_sra(n) = (-G_sra(n-1)/tau_sra)*dt + G_sra(n-1);
    if V(n) > V_th
        V(n) = V_reset;
        G_sra(n) = G_sra(n) + DeltaG_sra;
    end
    
end

subplot(3,1,1);
plot(t,I_app);
xlabel('Time/(s)');
ylabel('Current')

subplot(3,1,2);
plot(t,V);
xlabel('Time/(s)');
ylabel('Membrane Potential');

subplot(3,1,3);
plot(t,G_sra);
xlabel('Time/(s)');
ylabel('Adaption Conductance');


