% Aquesta funci� calcula el par�ntesi de Lie [f,g] de dos camps vectorials
% simb�lics f i g.
function h = Lie_p(f,g)
h = jacobian(g)*f - jacobian(f)*g;
end