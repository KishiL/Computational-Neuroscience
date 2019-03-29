function [sta,tcorr] = STA_spatial(Iapp, spikes, dt, tminus, tplus)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
% outputArg1 = inputArg1;
% tcorr = inputArg2;

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
    if spiketimes(n)-nminus > 0 && spiketimes(n)+nplus <= length(Iapp)
        sta = sta + Iapp(:,spiketimes(n)-nminus:spiketimes(n)+nplus);
    end
end

sta = sta/length(spiketimes);



end

