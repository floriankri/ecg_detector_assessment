function y = sortfilt1(x, n, p)

N = length(x);

if p>100
    p=100;
elseif p<0
    p=0;
end

if mod(n,2)==0 %even
    N1 = (n/2)-1;
    N2 = (n/2);
else
    N1 = (n-1)/2;
    N2 = (n-1)/2;
end

y = zeros(size(x));

USE_MEX = 1;
if USE_MEX
	if p == 0 %min filter
		y = minmaxfilter(x, n, 0);
	elseif p == 100 %max filter
		y = minmaxfilter(x, n, 1);
	elseif p == 50
		y = medianfilter(x, n);
	else
		for i=1:N
			A = max(1, i-N1);
			B = min(N, i+N2);
			P = 1 + round((p/100)*(B-A));
			Z = sort(x(A:B));
			y(i) = Z(P);
		end
	end
else
	for i=1:N
		A = max(1, i-N1);
		B = min(N, i+N2);
		P = 1 + round((p/100)*(B-A));
		Z = sort(x(A:B));
		y(i) = Z(P);
	end
end
