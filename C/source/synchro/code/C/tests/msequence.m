% MLS Sequence generator. Copyright Arthur C. Ludwig, March 2001

N = 5;

if N == 18; taps=[0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1]; end;

if N == 17; taps=[0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1]; end;

if N == 16; taps=[0 0 0 1 0 0 0 0 0 0 0 0 1 0 1 1]; end;

if N == 15; taps=[0 0 0 0 0 0 0 0 0 0 0 0 0 1 1]; end;

if N == 14; taps=[0 0 0 1 0 0 0 1 0 0 0 0 1 1]; end;

if N == 13; taps=[0 0 0 0 0 0 0 0 1 1 0 1 1]; end;

if N == 12; taps=[0 0 0 0 0 1 0 1 0 0 1 1]; end;

if N == 11; taps=[0 0 0 0 0 0 0 0 1 0 1]; end;

if N == 10; taps=[0 0 0 0 0 0 1 0 0 1]; end;

if N == 9; taps=[0 0 0 0 1 0 0 0 1]; end;

if N == 8; taps=[0 0 0 1 1 1 0 1]; end;

if N == 7; taps=[0 0 0 1 0 0 1]; end;

if N == 6; taps=[0 0 0 0 1 1]; end;

if N == 5; taps=[0 0 1 0 1]; end;

if N == 4; taps=[0 0 1 1]; end;

if N == 3; taps=[0 1 1]; end;

M = 2^N-1;

m = ones(1,N);
%disp(m);

regout = zeros(1,M);

for ind = 1:M

buf = mod(sum(taps.*m),2);
%disp(buf);

m(2:N) = m(1:N-1);

m(1)=buf;

%disp(m);

regout(ind) = m(N);

end

regout;