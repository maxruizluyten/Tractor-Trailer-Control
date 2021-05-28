% Aquesta funció calcula el polnomi interpolador de grau mínim que satisfà
% que la seva derivada i-èssima pren valors Z0(i) a temps T0 i ZF(i) a
% temps TF. La funció retorna un vector amb els coeficients del polinomi on
% el primer element és el terme de grau màxim i l'últim element és el
% temrme independent.
function coef = interpolate(Z0,ZF,T0,TF)
n = length(Z0);
syms t;

% Generem un vector simbolic p = [t, t, ..., t, 1]:
p = [1];
for i = 2:2*n
    p = [t, p];
end

% Usem el producte acomulat per igualar p = [t^{2n-1}, t^{2n-2}, ..., t, 1]
p = transpose(cumprod(p,'reverse'));
g(t) = transpose(p);

% Definim la matriu A del sistema d'equacions lineal que ens permet
% determinar els coeficients del polinomi interpolador:
A = [double(g(T0)); double(g(TF))];
for i = 2:n
    p(:,i) = jacobian(p(:,i-1));
    g(t) = transpose(p(:,i));
A = [A; double(g(T0)); double(g(TF))];
end

% Definim el vector b del sistema d'equacions lineals Ax = b:
b = [];
for i = 1:n
   b = [b; Z0(i); ZF(i)]; 
end

% Mostrem per pantalla el condicionament de A per assegurar que podem
% trobar solució al sistema:
disp(['El condicionament de la matriu és: ',num2str(cond(A))])

% Els coeficients del polinomi interpolador son la solució de Ax = b:
coef = mldivide(A,b);
end