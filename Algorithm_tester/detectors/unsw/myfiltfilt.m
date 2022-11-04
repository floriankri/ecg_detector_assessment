function y = myfiltfilt(b,a,x)

if length(x)<= 3*max([length(b)-1,length(a)-1])
    y=zeros(size(x));
else
    y=filtfilt(b,a,x);
end

return