% Aquesta funci� calcula la derivada de Lie d'una funci� simb�lica h
% respecte un camp vectorial simb�lic f.
function g = Lie_d(f,h)
g = jacobian(h)*f;
end