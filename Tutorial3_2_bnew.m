clear 
%Part 2a
rate = 20;
dt = 0.1E-3; % binsize
tmax = 100; %max time
tvec = (0:dt:tmax); %time vector

spikes = rand(1,length(tvec)-1)<rate*dt;

spike_times = find(spikes)*dt;
ISI = diff(spike_times);

figure(1);
histogram(ISI,25);
xlabel('ISI(sec)')
ylabel('frequency')

spikenew = expandspikebin(spikes, dt, 100E-3);
Fano = var(spikenew)./mean(spikenew);
cv = std(ISI)/mean(ISI);

k = 0;
for winsize = 10E-3:10E-3:1000E-3
    k = k+1;
    spikenew = expandspikebin(spikes, dt, winsize);
    Fano2(k) = var(spikenew)./mean(spikenew);
    window(k) = winsize;
end

figure(2);
plot(window,Fano2);
xlabel('window size');
ylabel('Fano Factor')




