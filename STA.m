function [sta, tcorr] = STA(Iapp, spikes, dt, tminus, tplus)
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
if(~exist('tminus'))
    tminus = 75E-3;
end

if(~exist('tplus'))
    tplus = 25E-3;
end

nminus = tminus/dt;
nplus = tplus/dt;
tcorr = (-nminus*dt):dt:(nplus*dt);
sta = zeros(size(tcorr));
spiketimes = find(spikes);

for n = 1:length(spiketimes)
    if spiketimes(n)-nminus > 0 && spiketimes(n)+nplus < length(Iapp)
        for k=1:length(sta)
            sta(k) = sta(k)+Iapp(spiketimes(n)-nminus+k-1);
        end

    end
end

sta = sta/sum(spiketimes);
