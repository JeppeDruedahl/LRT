function h = hermite(n,x)
% compute the Hermite polynomials.

x = x/sqrt(2);

% check n
if( n<0 ), error('The order of Hermite polynomial must be greater than or equal to 0.'); end

% again check n is an integer
if( 0~=n-fix(n) ), error('The order of Hermite polynomial must be an integer.'); end

% call the hermite recursive function.
h = hermite_rec(n);

% evaluate the hermite polynomial function, given x
if( nargin==2 )
    y = h(end) * ones(size(x));
    p = 1;
    for i=length(h)-1:-1:1
        y = y + h(i) * x.^p;
        p = p+1;
    end
    
    % restore the shape of y, the same as x
    h = reshape(y,size(x));
    h = 2^(-n/2)*h; 
end


function h = hermite_rec(n)
% This is the reccurence construction of a Hermite polynomial, i.e.:
%   H0(x) = 1
%   H1(x) = 2x
%   H[n+1](x) = 2x Hn(x) - 2n H[n-1](x)

if( 0==n ), h = 1;
elseif( 1==n ), h = [2 0];
else
    
    h1 = zeros(1,n+1);
    h1(1:n) = 2*hermite_rec(n-1);
    
    h2 = zeros(1,n+1);
    h2(3:end) = 2*(n-1)*hermite_rec(n-2);
    
    h = h1 - h2;
    
end