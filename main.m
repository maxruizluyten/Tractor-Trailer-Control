%% Definició del sistema
clear all
close all

% Dimensió del sistema prolongat:
n = 10;

% Donem valor a les constants de longitud del sistema:
d1 = 1;
d2 = 1;
d3 = 1;

% Definim les variables d'estat:
syms x;
syms y;
syms theta;
syms psi;
syms delta1;
syms delta2;
syms u10;
syms u11;
syms u12;
syms u2;

% Definim el vector q que conté totes les variables del sistema:
q = [x,y,theta,psi,delta1,delta2,u10,u11,u12,u2];


% Definim la funció auxiliar S que apareix a la quarta component de la 
% funció f, i les funcions auxiliars p barra i q barra que apareixen a
% l'expresió de S:
bar_p(q) = d2*sin(psi)*sin(delta1) + d1*cos(psi)*cos(delta1);
bar_q(q) = d2*cos(psi)*sin(delta1) - d1*sin(psi)*cos(delta1);
S(q) =  tan(delta1)/d1 - (bar_p*tan(delta2) - bar_q)/(d1*d3*cos(delta1));

% Definim les funcions que caracteritzen la dinàmica del nostre sistema:
f(q) = [cos(theta)*u10; sin(theta)*u10; tan(delta1)*u10/d1; S*u10; u2;
    0; u11; u12; 0; 0];
g1(q) = [x-x;0;0;0;0;0;0;0;1;0]; % escrivim x-x en comptes d'un 0 perquè
g2(q) = [x-x;0;0;0;0;0;0;0;0;1]; % el matlab detecti aquestes funcions com
g3(q) = [x-x;0;0;0;0;1;0;0;0;0]; % a funcions simbòliques de q.

%% Definició de condicions inicials i finals

% Fixem unes condicions inicials per al nostre sistema:
T0 = 0;

x_0 = 0;
y_0 = 0;
theta_0 = pi/4;
psi_0 = 0;
delta1_0 = 0;
delta2_0 = 0;
u10_0 = 1;
u11_0 = 0;
u12_0 = 0;
u2_0 = 0;

q0 = [x_0, y_0, theta_0, psi_0, delta1_0, delta2_0, u10_0, u11_0,...
    u12_0, u2_0];

% Fixem unes condicions finsals a les quals volen que arribi el nostre
% sistema:
TF =20;

x_F = 20;
y_F = 20;
theta_F = pi/4;
psi_F = pi/10;
delta1_F = 0;
delta2_F = 0;
u10_F = 1;
u11_F = 0;
u12_F = 0;
u2_F = 0;

qF = [x_F, y_F, theta_F, psi_F, delta1_F, delta2_F, u10_F, u11_F,...
    u12_F, u2_F];

% Definim les tres funcions h1, h2 i h3 per trobar el canvi de variables a
% la forma canònica de Brunovsky:
h1(q) = x;
h2(q) = y;
h3(q) = psi;

% Usem l'algorisme vist a teoria per calcular les noves variables:
z1(q) = h1;
z2(q) = Lie_d(f,z1);
z3(q) = Lie_d(f,z2);
z4(q) = Lie_d(f,z3);
z5(q) = h2;
z6(q) = Lie_d(f,z5);
z7(q) = Lie_d(f,z6);
z8(q) = Lie_d(f,z7);
z9(q) = h3;
z10(q) = Lie_d(f,z9);

% Calculem les corresponents condicions inicials amb les variables z:
z1_0 = double(z1(q0(1), q0(2), q0(3), q0(4), q0(5), q0(6), q0(7),...
    q0(8), q0(9), q0(10)));
z2_0 = double(z2(q0(1), q0(2), q0(3), q0(4), q0(5), q0(6), q0(7),...
    q0(8), q0(9), q0(10)));
z3_0 = double(z3(q0(1), q0(2), q0(3), q0(4), q0(5), q0(6), q0(7),...
    q0(8), q0(9), q0(10)));
z4_0 = double(z4(q0(1), q0(2), q0(3), q0(4), q0(5), q0(6), q0(7),...
    q0(8), q0(9), q0(10)));
z5_0 = double(z5(q0(1), q0(2), q0(3), q0(4), q0(5), q0(6), q0(7),...
    q0(8), q0(9), q0(10)));
z6_0 = double(z6(q0(1), q0(2), q0(3), q0(4), q0(5), q0(6), q0(7),...
    q0(8), q0(9), q0(10)));
z7_0 = double(z7(q0(1), q0(2), q0(3), q0(4), q0(5), q0(6), q0(7),...
    q0(8), q0(9), q0(10)));
z8_0 = double(z8(q0(1), q0(2), q0(3), q0(4), q0(5), q0(6), q0(7),...
    q0(8), q0(9), q0(10)));
z9_0 = double(z9(q0(1), q0(2), q0(3), q0(4), q0(5), q0(6), q0(7),...
    q0(8), q0(9), q0(10)));
z10_0 = double(z10(q0(1), q0(2), q0(3), q0(4), q0(5), q0(6), q0(7),...
    q0(8), q0(9), q0(10)));

% Calculem les corresponents condicions inicials amb les variables z:
z1_F = double(z1(qF(1), qF(2), qF(3), qF(4), qF(5), qF(6), qF(7),...
    qF(8), qF(9), qF(10)));
z2_F = double(z2(qF(1), qF(2), qF(3), qF(4), qF(5), qF(6), qF(7),...
    qF(8), qF(9), qF(10)));
z3_F = double(z3(qF(1), qF(2), qF(3), qF(4), qF(5), qF(6), qF(7),...
    qF(8), qF(9), qF(10)));
z4_F = double(z4(qF(1), qF(2), qF(3), qF(4), qF(5), qF(6), qF(7),...
    qF(8), qF(9), qF(10)));
z5_F = double(z5(qF(1), qF(2), qF(3), qF(4), qF(5), qF(6), qF(7),...
    qF(8), qF(9), qF(10)));
z6_F = double(z6(qF(1), qF(2), qF(3), qF(4), qF(5), qF(6), qF(7),...
    qF(8), qF(9), qF(10)));
z7_F = double(z7(qF(1), qF(2), qF(3), qF(4), qF(5), qF(6), qF(7),...
    qF(8), qF(9), qF(10)));
z8_F = double(z8(qF(1), qF(2), qF(3), qF(4), qF(5), qF(6), qF(7),...
    qF(8), qF(9), qF(10)));
z9_F = double(z9(qF(1), qF(2), qF(3), qF(4), qF(5), qF(6), qF(7),...
    qF(8), qF(9), qF(10)));
z10_F = double(z10(qF(1), qF(2), qF(3), qF(4), qF(5), qF(6), qF(7),...
    qF(8), qF(9), qF(10)));

%% Determinació del sistema i els controls desitjats en les variables z

% Per calcular els controls desitjats w (en les variables z), usarem
% interpolació polinòmica. Usarem la notació z_s per designar les funcions
% del temps z que son solució del problema amb les condicions inicials i
% finals fixades prèviament.

% Definim una variable simbòlica t per a designar el temps: 
syms t;

% Calculem el primer control w1 i les primeres quatre variables z1_s, z2_s,
% z3_s i z4_s:
Z0 = [z1_0, z2_0, z3_0, z4_0];
ZF = [z1_F, z2_F, z3_F, z4_F];
coef = interpolate(Z0,ZF,T0,TF);
p(t) = t-t; % Escribim t-t en comptes de 0 perquè el Matlab detecti p com a
            % polinomi simbòlic.
l = length(coef); % Coef guarda els coeficients del polinomi interpolador

% Generem el polinomi interpolador:
for i = 1:l
    p(t) = p + coef(i)*t^(l-i); 
end

% Calculem les variables solució i el control w1:
z1_s(t) = p;
z2_s(t) = jacobian(z1_s);
z3_s(t) = jacobian(z2_s);
z4_s(t) = jacobian(z3_s);
w1(t) = jacobian(z4_s);

% Calculem el segon control w2 i les primeres quatre variables z5_s, z6_s,
% z7_s i z8_s de manera anàloga:
Z0 = [z5_0, z6_0, z7_0, z8_0];
ZF = [z5_F, z6_F, z7_F, z8_F];
coef = interpolate(Z0,ZF,T0,TF);
p(t) = t-t; % Sobreescrivim el polinomi p
l = length(coef);
for i = 1:l
    p(t) = p + coef(i)*t^(l-i); 
end
z5_s(t) = p;
z6_s(t) = jacobian(z5_s);
z7_s(t) = jacobian(z6_s);
z8_s(t) = jacobian(z7_s);
w2(t) = jacobian(z8_s);

% Calculem el tercer control w3 i les primeres quatre variables z5_s, z6_s
% de manera anàloga:
Z0 = [z9_0, z10_0];
ZF = [z9_F, z10_F];
coef = interpolate(Z0,ZF,T0,TF);
p(t) = t-t; % Sobreescrivim el polinomi p
l = length(coef);
for i = 1:l
    p(t) = p + coef(i)*t^(l-i); 
end
z9_s(t) = p;
z10_s(t) = jacobian(z9_s);
w3(t) = jacobian(z10_s);

%% Càlcul de la solució del sistema. Mètode 1:
% En aquesta secció, desfarem la realimentació que
% hem usat per calcular el control que ens porta les varibles del
% sistema del punt inicial al punt final desitjats. Posteriorment,
% resoldrem la equació diferencial del sistema numèricament per trobar les
% variables d'estat.

% w és el control en les varibles z. v és el control en les varibales q. La
% realimentació està definida com w = a + b*v. Definim a i b i coloquem w
% com a vector:
w(t) = [w1; w2; w3];
a(q) = [Lie_d(f,z4); Lie_d(f,z8); Lie_d(f,z10)];
b(q) = [Lie_d(g1,z4), Lie_d(g2,z4), Lie_d(g3,z4); Lie_d(g1,z8),... 
        Lie_d(g2,z8), Lie_d(g3,z8); Lie_d(g1,z10), Lie_d(g2,z10),...
        Lie_d(g3,z10)];
    
% Escrivim g com a vector:
g(q) = [g1, g2, g3];

% Aleshores el sistema ve definit per: 
% dx/dt = f + g*v = f + g*(w-a)*inv(b)
% Creem una funció numèrica F = f + g*(w-a)*inv(b) per resoldre la edo:
F = @(s,X) double(f(X(1),X(2),X(3),X(4),X(5),X(6),X(7),X(8),X(9),X(10))...
    -g(X(1),X(2),X(3),X(4),X(5),X(6),X(7),X(8),X(9),X(10))*...
    inv(b(X(1),X(2),X(3),X(4),X(5),X(6),X(7),X(8),X(9),X(10)))*...
    a(X(1),X(2),X(3),X(4),X(5),X(6),X(7),X(8),X(9),X(10))+...
    g(X(1),X(2),X(3),X(4),X(5),X(6),X(7),X(8),X(9),X(10))*...
    inv(b(X(1),X(2),X(3),X(4),X(5),X(6),X(7),X(8),X(9),X(10)))*...
    w(s));

% Resolem la equació deiferencial amb la funció ode45:
[s,Q_metode1] = ode45(F,[T0,TF],q0);

figure(1)
set(gcf,'color','w')
subplot(2,5,1)
plot(s,Q_metode1(:,1))
title('Evolution of $x$ in time','Interpreter','Latex')
xlabel('$t$','Interpreter', 'Latex') 
ylabel('$x(t)$','Interpreter', 'Latex')

subplot(2,5,2)
plot(s,Q_metode1(:,2))
title('Evolution of $y$ in time','Interpreter','Latex')
xlabel('$t$','Interpreter', 'Latex') 
ylabel('$y(t)$','Interpreter', 'Latex')

subplot(2,5,3)
plot(s,Q_metode1(:,3))
title('Evolution of $\theta$ in time','Interpreter','Latex')
xlabel('$t$','Interpreter', 'Latex') 
ylabel('$\theta(t)$','Interpreter', 'Latex')

subplot(2,5,4)
plot(s,Q_metode1(:,4))
title('Evolution of $\psi$ in time','Interpreter','Latex')
xlabel('$t$','Interpreter', 'Latex') 
ylabel('$\psi(t)$','Interpreter', 'Latex')

subplot(2,5,5)
plot(s,Q_metode1(:,5))
title('Evolution of $\delta_1$ in time','Interpreter','Latex')
xlabel('$t$','Interpreter', 'Latex') 
ylabel('$\delta_1(t)$','Interpreter', 'Latex')

subplot(2,5,6)
plot(s,Q_metode1(:,6))
title('Evolution of $\delta_2$ in time','Interpreter','Latex')
xlabel('$t$','Interpreter', 'Latex') 
ylabel('$\delta_2(t)$','Interpreter', 'Latex')

subplot(2,5,7)
plot(s,Q_metode1(:,7))
title('Evolution of $u_1^0$ in time','Interpreter','Latex')
xlabel('$t$','Interpreter', 'Latex') 
ylabel('$u_1^0(t)$','Interpreter', 'Latex')

subplot(2,5,8)
plot(s,Q_metode1(:,8))
title('Evolution of $u_1^1$ in time','Interpreter','Latex')
xlabel('$t$','Interpreter', 'Latex') 
ylabel('$u_1^1(t)$','Interpreter', 'Latex')

subplot(2,5,9)
plot(s,Q_metode1(:,9))
title('Evolution of $u_1^2$ in time','Interpreter','Latex')
xlabel('$t$','Interpreter', 'Latex') 
ylabel('$u_1^2(t)$','Interpreter', 'Latex')

subplot(2,5,10)
plot(s,Q_metode1(:,10))
title('Evolution of $u_2$ in time','Interpreter','Latex')
xlabel('$t$','Interpreter', 'Latex') 
ylabel('$u_2(t)$','Interpreter', 'Latex')

%% Càlcul de la solució del sistema. Mètode 2:
% En aquesta secció invertirem el canvi de variables per trobar les
% variables d'estat en funció del temps que resolen el problema de
% condicions inicials i finals. Posteriorment graficarem aquestes
% variables.

% Aquestes són les variables solució. Les expressions de p barra i q barra
% en funció del temps s'inclouen per poder definir delta2.
x_s(t) = z1_s;
y_s(t) = z5_s;
theta_s(t) = atan(z6_s/z2_s);
psi_s(t) = z9_s;
delta1_s(t) = atan(d1*(-z3_s*sin(theta_s)+z7_s*cos(theta_s))/...
    (z2_s^2+z6_s^2));
u10_s(t) = sqrt(z2_s^2 + z6_s^2);
bar_p_s(t) = d2*sin(psi_s)*sin(delta1_s) + d1*cos(psi_s)*cos(delta1_s);
bar_q_s(t) = d2*cos(psi_s)*sin(delta1_s) - d1*sin(psi_s)*cos(delta1_s);
delta2_s(t) = atan(-d1*d3*cos(delta1_s)/bar_p_s*(z10_s/u10_s - ...
    (d1*(-z3_s*sin(theta_s)+z7_s*cos(theta_s))/(z2_s^2+z6_s^2))...
    /d1 - bar_q_s/(d1*d3*cos(delta1_s))));
u11_s(t) = jacobian(u10_s);
u12_s(t) = jacobian(u11_s);
u2_s(t) = jacobian(delta1_s);

q_s(t) = [x_s(t), y_s(t), theta_s(t), psi_s(t), delta1_s(t),...
    delta2_s(t), u10_s(t), u11_s(t), u12_s(t), u2_s(t)];

% Grafiquem les funcions obtingudes (de fet, només en grafiquem 4 de les 10
% funcions):
s = linspace(T0,TF,100);
Q_metode2 = zeros(100,10);
for j = 1:length(s)
    Q_metode2(j,:) = double(q_s(s(j)));
end
figure(2)
set(gcf,'color','w')
subplot(2,2,1)
plot(s,Q_metode2(:,1))
title('Evolution of $x$ in time','Interpreter','Latex')
xlabel('$t$','Interpreter', 'Latex') 
ylabel('$x(t)$','Interpreter', 'Latex')

subplot(2,2,2)
plot(s,Q_metode2(:,4))
title('Evolution of $\psi$ in time','Interpreter', 'Latex')
xlabel('$t$','Interpreter', 'Latex') 
ylabel('$\psi(t)$','Interpreter', 'Latex')

subplot(2,2,3)
plot(s,Q_metode2(:,7))
title('Evolution of $u_1^0$ in time','Interpreter', 'Latex')
xlabel('$t$','Interpreter', 'Latex') 
ylabel('$u_1^0(t)$','Interpreter', 'Latex')

subplot(2,2,4)
plot(s,Q_metode2(:,9))
title('Evolution of $u_1^2$ in time','Interpreter', 'Latex')
xlabel('$t$','Interpreter', 'Latex') 
ylabel('$u_1^2(t)$','Interpreter', 'Latex')

% Grafiquem també les controls que usem per solucionar el sistema:
v1_s(t) = jacobian(u12_s);
v2_s(t) = jacobian(u2_s);
v3_s(t) = jacobian(delta2_s);

for j = 1:length(s)
    control_1(j) = double(v1_s(s(j)));
end
for j = 1:length(s)
    control_2(j) = double(v2_s(s(j)));
end
for j = 1:length(s)
    control_3(j) = double(v3_s(s(j)));
end

figure(3)
set(gcf,'color','w')
subplot(3,1,1)
plot(s,control_1)
title('Evolution of $v_1$ in time','Interpreter','Latex')
xlabel('$t$','Interpreter', 'Latex') 
ylabel('$v_1(t)$','Interpreter', 'Latex')

subplot(3,1,2)
plot(s,control_2)
title('Evolution of $v_2$ in time','Interpreter', 'Latex')
xlabel('$t$','Interpreter', 'Latex') 
ylabel('$v_2(t)$','Interpreter', 'Latex')

subplot(3,1,3)
plot(s,control_3)
title('Evolution of $v_3$ in time','Interpreter', 'Latex')
xlabel('$t$','Interpreter', 'Latex') 
ylabel('$v_3(t)$','Interpreter', 'Latex')