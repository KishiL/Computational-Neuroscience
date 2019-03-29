function [new_vector] = expandbin(old_vector, old_dt, new_dt)
len_oldvec = length(old_vector)-1;
ratio = new_dt/old_dt;
len_newvec = ceil(len_oldvec/ratio);
new_vector = zeros(1,len_newvec);

for k= 1:length(new_vector)
    a = mean(old_vector((k-1)*ratio+1:k*ratio));
    if a > 0.01
        new_vector(k) = 1;
    else
        new_vector(k) = a;
    end
end

end

