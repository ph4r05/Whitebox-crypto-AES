
% just simple number of bases

% parameter defining matrix size
p = 8;

s=1;
for i=0:(p-1)
	s*=2^p-2^i;
end

total = 2^(p*p);
probab = s/total;
printf ("Invertible: %ld;\nTotal:      %ld\nProbability that we will succeed finding invertible matrix is: %f\n\n", s, total, probab)


% now compute expectation up to the K iterations
K=1000
EX=0
for i=K:-1:1
	curProbab=((1-probab)^(i-1)) * probab;
	EX+=i*curProbab;
end
printf ("Expectation for %d rounds: %f\n", K, EX)


