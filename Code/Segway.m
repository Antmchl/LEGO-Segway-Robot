% Work of MICHEL Antoine and Mounier Timothée
% November 2023
%% Physical constants
g = 9.81;                       % gravity acceleration [m/sec^2]
%% Physical parameters
m =    23*10^-3  ;						% wheel weight [kg]
R = 0.027;						% wheel radius [m]
Jw = m * R^2 / 2;				% wheel inertia moment [kgm^2]
M =   (617-2*23)*10^-3 ;                        % body weight [kg]
W = 0.105;						% body width [m]
D = 0.1;						% body depth [m]
h = 0.21;       				% body height [m]
L = 0.12;						% distance of the center of mass from the wheel axle [m]
Jpsi = M * L^2 / 3;				% body pitch inertia moment [kgm^2]
Jphi = M * (W^2 + D^2) / 12;	% body yaw inertia moment [kgm^2]
fm = 0.0022;					% friction coefficient between body & DC motor
fw = 0;           				% friction coefficient between wheel & floor
%% Motors parameters
Jm = 1e-5;						% DC motor inertia moment [kgm^2]
Rm = 6.69;						% DC motor resistance [Om]
Kb = 0.468;						% DC motor back EMF constant [Vsec/rad]
Kt = 0.317;						% DC motor torque constant [Nm/A]
K_PWM = 8.087;                  % Volts to PWM value coefficient [1/V]


%% Constantes definition
alpha =  Kt / Rm;
beta =  Kt * Kb / Rm + fm;
E_11 = (2 * m + M) * R^2 + 2 * Jw + 2 * Jm;
E_12 = M * L * R - 2  * Jm;
E_22 = M * L^2 + Jpsi + 2 * Jm;


%% Part 1: Model

% définition des matrices selon la préparation, on transpose le vecteur pour
% H
E=[E_11 E_12;
    E_12 E_22];

F=[ 2*beta -2*beta;
    -2*beta 2*beta];

G= [0 0;
    0 -M*g*L];

H= [2*alpha -2*alpha]';

% création d'un vecteur O pour la création de la matrice B1
O = [0;
    0];

A1= [ 0 0 1 0;
    0 0 0 1;
    -(inv(E))*G -(inv(E))*F];

B1=[0 ;
    0 ;
    (inv(E))*H ];

C=eye(4);

C1 =ss(A1,B1,C,0);
s1.StateName = {'theta', 'psi', 'theta_dot', 'psi_dot'};
s1.InputName = {'U'};


% 3) Les états d'équilibres correspondent aux états où
% Psi_dot, theta_dot et Psi sont nuls
% On voit au près de la matrice A que l'on a un deg de liberté en theta
% Car on un une colonne de 0 (chute du rang)

% 4) a) On place le segway devant nous et on le lache 
% suspens ... il tombe

% b) on calcul avec pzmap on trouve les valeurs de pôles contenues
% dans p [0 -370.06 6.83 -6.50] les pôles ne sont pas tous a partie 
% réelle strict négative donc le système est instable

[p,z]= pzmap(C1);



%% For discrete control and simulation
Ts = 0.004;                     % Control system sample time
Psi0 = deg2rad(10);             % Initial value to disturb the system

% bloqueur d'ordre Zéro car valeur constante
sysd=c2d(C1,Ts,'zoh');
Phi=sysd.a;
Gamma=sysd.b;
Cd=sysd.c;

% vérifiaction que phi(ts)=exp(A1*Ts) et on obtient bien la même chose;

Phi_Ts = expm(A1.*Ts);


C2 =ss(Phi,Gamma,Cd,0);                        

% En effet quand on réalise les calculs suivants on peut confirmer qu'il
% existe un lien entre pole continu et discret
% la limite de stabilité correspond au cercle unité dans le complexe car 
% il s'agit d'une expo complexe.                                                                                                                                                                                                           
% les poles stables sont à l'intérieur strictement.
[pd,zd]=pzmap(C2);
zp=exp(p(:).*Ts);


% Q4)
% on obtient la matrice de commandabilité en discret avec phi et gamma
% et on trouve bien une matrice de rang 4 donc la matrice est commandable
% si elle n'était pas commandable, on ne pourrait pas stabiliser le système
Command=[Gamma Phi*Gamma (Phi^2)*Gamma (Phi^3)*Gamma];
rank(Command);

%% Part 2

% Q1)
% Xeq va donc tendre vers une valeur nulle car on a une suite géo dont
% la raiqon en valeur abs est inf stricte à 1
% (qui n'est pas une abération) 
% On va choisir des pôles positif pour éviter d'avoir une suite alternée
% qui abimerait dans le temps la motorisation


% Q2) 
% le pôle dominant correspond aux constantes de temps les plus élevées ie le
% plus lent,il rend compte de la dynamique du système
pd;
Td=1./pd;

% on a comme pole stables : 0.23 et 0.97
% et les constantes de temps associées sont 4.39 et 1.03
% un pole lent se situe au plus proche du cercle unité

% Q3)
x=0.9845;
P=[pd(2) pd(4) x x];
K=acker(Phi,Gamma,P);
% quand on prend une valeur proche de 0.97 on se trouve à une tension
% supérieure à 10V mais quand on tend vers 1 la tension est inférieure de
% 10V

% analyse des courbes : stabilité des états, on voit bien l'avancée du
% segway et le recul
% de plus pour la tension on voit bien le recul avec la tension qui passe ne
% négatif

 

%% partie 6.1
% gain K_Pwm
% round justifie la résolution des capteurs



% Part 6.2

% Q2)
f_oub=0.999;
% Y(z)=f_oub*(z^-1)*Y(z) + (1-f_oub)*U(z)

numerator = [(1-f_oub) 0];
denominator = [1 -f_oub];
sysf = tf(numerator,denominator,Ts);

bode(sysf)

% wc = 0.25; %rad/s
% Oui c'est le filtre adéquat pour détecter notre biais car
% il va être visible sur notre BODE pour une freq nulle

% Part 7.1

% la variable du vecteur d'état qu 'il est possible d'imposer est psi

% pour les saccades il s'agissait d'une erreur avec K_PWM

Q = inv(eye(4) + Gamma*K - Phi);
G = [1 0 0 0];
N = inv(G*Q*Gamma);

% Q3)
% on arrive bien à le faire avancer et reculer avec un seuil  
% défini en cm.
