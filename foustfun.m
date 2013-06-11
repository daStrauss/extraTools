N = 512;
Nsamp = ceil(1.5*5*log(N));

% Construct half the FFT matrix
W = zeros(N/2,N);
w = exp(-2*pi*i/N);
for ii=1:N/2
 for jj=1:N
   W(ii,jj) = w^((ii-1)*(jj-1));
 end
end

% Sample f
x = linspace(-1,1,N).';
% frequencies
freq = [10 4 2 13 9];
amps = [.5 .2 .9 3 1];
fall = sum( (ones(N,1)*amps) .* sin(2*pi*x*freq*(N-1)/N),2 );

% Sample the signal at random places
isamp = floor(rand(Nsamp-1,1)*N)+1;
fsamp = fall(isamp);

% Reconstruct
cvx_begin
 variable g(N)
 minimize norm(W*g,1)
 subject to
   g(isamp) == fsamp
cvx_end

figure(100)
clf
subplot(311)
plot(fall);
title('original with sample locations')
hold on
plot(isamp,fall(isamp),'x');
subplot(312)
plot(g);
title('reconstruction');

subplot(313)
plot(1:256, abs(W*fall), 1:256, abs(W*g), 'o')
