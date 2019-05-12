clear 
%Part 2b
rate = 20;
dt = 0.1E-3; % binsize
tmax = 10; %max time
tvec = (0:dt:tmax); %time vector
spikes = zeros(1000,length(tvec));

for i = 1:1000
    spikes(i,:) = rand(1,length(tvec))<rate*dt;
end

csspikes = cumsum(spikes,2);
Fano = var(csspikes)./mean(csspikes);
ax1 = subplot(2,1,1);
xlabel('time (sec)')
ylabel('Fano Factor')
plot(ax1,tvec,Fano);

% alternative Fano Factor computation
seg_len = round(0.2/dt);
for trial = 1:1000
    for m = 1:round(tmax/0.2)
        spike_cnt(trial,m) = sum(spikes(trial,(m-1)*seg_len+1:m*seg_len));
    end
end
Fano2 = var(spike_cnt)./mean(spike_cnt);
ax2 = subplot(2,1,2)
xlabel('time (sec)')
ylabel('Fano Factor')
plot(ax2,[0:0.2:tmax-0.2],Fano2)

axis([ax1,ax2],[0 10 0.9 1.15])
title(ax1,'Method 1')
title(ax2,'Alternative Method')

%inhomogeneous Poisson process
ratet = 25 + 20*sin(2*pi*tvec);
for k = 1:1000
    ihspikes(k,:) = rand(size(tvec))< ratet*dt;
end

csspikes_new = cumsum(ihspikes,2);
Fano3 = var(csspikes_new)./mean(csspikes_new);
figure(3);
plot(tvec,Fano3);
xlabel('time (sec)')
ylabel('Fano Factor')

%alternative Fano Factor computation
seg_len = round(0.2/dt);
for trial = 1:1000
    for m = 1:round(tmax/0.2)
        spike_cnt(trial,m) = sum(ihspikes(trial,(m-1)*seg_len+1:m*seg_len));
    end
end
Fano4 = var(spike_cnt)./mean(spike_cnt);
figure(4)
plot([0:0.2:tmax-0.2],Fano4)
xlabel('time (sec)')
ylabel('Fano Factor')



