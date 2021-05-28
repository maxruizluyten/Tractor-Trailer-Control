% Aquesta funció calcula la derivada de Lie d'una funció simbòlica h
% respecte un camp vectorial simbòlic f.
function g = Lie_d(f,h)
g = jacobian(h)*f;
end