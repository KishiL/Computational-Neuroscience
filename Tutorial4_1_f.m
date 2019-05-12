clear
%
G_L = 30E-9;
G_Na_max = 12E-6;
G_K_max = 3.6E-6;
E_Na = 45E-3;
E_K = -82E-3;
E_L = -0.065;
E_L_rev = -60E-3;
C_m = 100E-12;

dt = 5E-7; t_max = 0.35;
t = 0:dt:t_max;

% Get the I_app required: baseline 0.6nA, 
I_0 = 0.7E-9;
I_app = ones(size(t))*I_0;
I_app(floor(0.1/dt): floor((0.1+0.005)/dt)) = 1E-9;

V_m = zeros(size(t));
V_m(1) = E_L;
m = zeros(size(t));
h = zeros(size(t));
n = zeros(size(t));

m(1) = 0;
h(1) = 0;
n(1) = 0;

for k = 2:length(t)-1
    
    alpha_m = (1E5*(-V_m(k-1)-0.045))/(exp(100*(-V_m(k-1) - 0.045))-1);
    beta_m = 4*(1E3)*exp((-V_m(k-1) - 0.070)/0.018);
    alpha_h = 70*exp(50*(-V_m(k-1)-0.070));
    beta_h = (1E3)/(1+exp(100*(-V_m(k-1) - 0.040)));
    alpha_n = (1E4)*(-V_m(k-1) - 0.060)/(exp(100*(-V_m(k-1)-0.060))-1);
    beta_n = 125*exp((-V_m(k-1) - 0.070)/0.08);
    
    %Model
    V_m(k) = (1/C_m)*(G_L*(E_L_rev - V_m(k-1)) + G_Na_max * (m(k-1))^3*h(k-1)*(E_Na-V_m(k-1)) + G_K_max*n(k-1)^4*(E_K - V_m(k-1)) + I_app(k-1))*dt + V_m(k-1);
    %    V_m(k) = (1/C_m)*(G_L*(E_L - V_m(k-1)) + G_Na_max * (m_std)^3*h_std*(E_Na-V_m(k-1)) + G_K_max * (n_std)^4 * (E_K - V_m(k-1)) + I_app(k-1))*dt + V_m(k-1);
    
    m(k) = (alpha_m*(1-m(k-1))-beta_m*m(k-1))*dt + m(k-1);
    h(k) = (alpha_h*(1-h(k-1))-beta_h*h(k-1))*dt + h(k-1);
    n(k) = (alpha_n*(1-n(k-1))-beta_n*n(k-1))*dt + n(k-1);
    
end

figure(1);
subplot(2,1,1);
plot(I_app);
ylabel('I_{app}/A')
subplot(2,1,2);
plot(V_m);
xlabel('t')
ylabel('Vm/V')



