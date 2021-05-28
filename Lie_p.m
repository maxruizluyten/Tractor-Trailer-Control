% Aquesta funció calcula el parèntesi de Lie [f,g] de dos camps vectorials
% simbòlics f i g.
function h = Lie_p(f,g)
h = jacobian(g)*f - jacobian(f)*g;
end