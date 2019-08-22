function overlap_nmi = ONMI(X, Y)

%%
% X: N x KX; Y: N x KY;
% Detecting the overlapping and hierarchical
% community structure in complex networks
%%

[~, KX] = size(X);
KY = size(Y, 2);

X = X';
Y = Y';

HX = HXi(X);
HY = HXi(Y);


HX_Y = 0;
for i = 1:KX
    tmp = HXiY(X(i, :), Y);
    if tmp == Inf
        tmp = HX(i);
    end
    tmp2 = HX(i);
    if tmp2 == 0
        tmp2 = 1;
    end
    HX_Y = HX_Y + tmp / tmp2;
end
HX_Y = HX_Y / KX;

HY_X = 0;
for i = 1:KY
    tmp = HXiY(Y(i, :), X);
    if tmp == Inf
        tmp = HY(i);
    end
    tmp2 = HY(i);
    if tmp2 == 0
        tmp2 = 1;
    end
    HY_X = HY_X + tmp / tmp2;
end
HY_X = HY_X / KY;

overlap_nmi = 1 - 0.5 * (HX_Y + HY_X);


end

function hx_i = HXi(X)
[kx, n] = size(X);
hx_i = zeros(1, kx);
for i = 1:kx
    xi = X(i, :);
    cnt = sum(xi);
    hx_i(i) = logEn(cnt, n) + logEn(n - cnt, n);
end
end


function ve = HXiY(Xi, Y)
ve = Inf;
m = size(Y, 1);
for j = 1:m
    Yj = Y(j, :);
    xi_yj = condition_entropy(Xi, Yj);
    ve = min(ve, xi_yj);
end
end




function x_con_y = condition_entropy(x, y)
n = length(x);
a = 0;
b = 0;
c = 0;
d = 0;
for i = 1:n
    if x(i) == 0 && y(i) == 0
        a = a + 1;
    elseif x(i) == 0 && y(i) == 1
        b = b + 1;
    elseif x(i) == 1 && y(i) == 0
        c = c + 1;
    else
        d = d + 1;
    end   
end

an = logEn(a, n);
bn = logEn(b, n);
cn = logEn(c, n);
dn = logEn(d, n);
bdn = logEn(b + d, n);
acn = logEn(a + c, n);

if an + dn > bn + cn
    x_con_y = an + bn + cn + dn - bdn - acn;
else
    x_con_y = Inf;
end
end


function wn = logEn(w, n)
if w <= 0
    wn = 0;
else
    p = w * 1.0 / n; 
    wn = -p * log2(p);
end
end