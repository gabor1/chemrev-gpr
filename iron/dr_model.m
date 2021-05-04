eVtoGPa = 160.21766208;

%%% GAP data

data_gap

Ek = 0.99473;
DEp_b = 49.45/1000;

%%% Dorn and Rajnak 1964 model

% Input

% Input Dorn/Rajnak Table 1
%R = 1.01;		% Gamma_c / Gamma_0
%alpha = 0.0010; 	% Peierls potential exponent
%El = 0.54201;		% Elastic energy

% Fit line tension / scaled plot
%R = 1.0069;
%alpha = -0.25;

% Fit line tension / total plot
R = 1.005638;
alpha = -0.5965;

% Easy model / Input alpha = 0
%R = 1.0071;
%alpha = 0;


% Numerical parameters
Nsteps = 10000;		% Number of integration steps for Un
Nsteps_k = 10000;	% Number of integration steps for Uk

% Primary quantities
b=2.4552;
a_b = 2*b/sqrt(3);
a = sqrt(2/3)*a_b;

DEp = DEp_b/b;

Gamma_0 = DEp/(R-1);

% Derived quantities
Gamma_c = Gamma_0*R;

if (abs(alpha)~=0)
tau_P = pi*(Gamma_c - Gamma_0)/(16*abs(alpha)*a*b)*(3+sqrt(1+8*alpha^2))*sqrt(8*alpha^2-2+2*sqrt(1+8*alpha^2));
else
tau_P = (pi*(Gamma_c-Gamma_0))/(a*b);
end

%K = [0.01 0.08030 0.15060 0.22090 0.29120];
K = [0.001, 0.01:0.01:0.99, 0.999];
y0 = [];

lagr = [];
Un = [];
Un_d = [];

% Check Peierls potential

Gamma_check = [];
x_step=-a/2:0.01:a/2;

for i=1:size(x_step,2)
	Gamma_check(i) = (Gamma_c+Gamma_0)/2+(Gamma_c-Gamma_0)/2*(alpha/4+cos(2*pi*x_step(i)/a)-alpha/4*cos(4*pi*x_step(i)/a));

end

for i=1:size(x_step,2)
	Gamma_check2(i) = (Gamma_c+Gamma_0)/2+(Gamma_c-Gamma_0)/2*(cos(4*pi*x_step(i)/a));

end

%Peierls_check = [];
%x_step=-a/2:0.01:a/2;

%for i=1:size(x_step,2)

%Peierls_check(i) = DEp/4*(3/2-2*cos(2*pi*x_step(i)/a)-(cos(2*pi*x_step(i)/a))^2);

%end


% Step 0 - Integrate numerically the kink self energy (reference energy)
Uk = 0;
Dy_k = a/Nsteps_k;

for j=1:Nsteps_k
	y01_k(j) = -a/2+(j-1)*Dy_k;
	y02_k(j) = -a/2+j*Dy_k;

	Gamma_y01_k(j) = (Gamma_c+Gamma_0)/2+(Gamma_c-Gamma_0)/2*(alpha/4+cos(2*pi*y01_k(j)/a)-alpha/4*cos(4*pi*y01_k(j)/a));
	Gamma_y02_k(j) = (Gamma_c+Gamma_0)/2+(Gamma_c-Gamma_0)/2*(alpha/4+cos(2*pi*y02_k(j)/a)-alpha/4*cos(4*pi*y02_k(j)/a));

	A1(j) = sqrt((Gamma_y01_k(j)/Gamma_0).^2-1);
	A2(j) = sqrt((Gamma_y02_k(j)/Gamma_0).^2-1);

	Uk = Uk + Gamma_0/2*Dy_k*(A1(j) + A2(j));
end

for i=1:size(K,2)

	% Step 1 - Calculate y0 (equilibrium position of straight dislocation / standard equilibrium)
	fun = @(y) K(i)*tau_P*b + pi/a*(Gamma_c-Gamma_0)*sin(2*pi.*y/a).*(1-alpha*cos(2*pi.*y/a));
	y0(i) = fzero(fun,-a/2);

	% Step 2 - Calculate lambda_c (critical position of the kink, s.t. dy/dx = / to avoid force concentrations)
	fun2 = @(lambda) (Gamma_c-Gamma_0)/2*(cos(2*pi.*lambda/a) - cos(2*pi*y0(i)/a)) - alpha/8*(Gamma_c - Gamma_0)*(cos(4*pi.*lambda/a)-cos(4*pi*y0(i)/a))-K(i)*tau_P*b*(lambda-y0(i));
	lambda_c(i) = fzero(fun2,0.01);

	% Step 3 - Integrate numerically the kink-pair nucleation energy
	Dy(i) = (lambda_c(i) - y0(i))/Nsteps;

	Gamma_y0(i) = (Gamma_c+Gamma_0)/2+(Gamma_c-Gamma_0)/2*(alpha/4+cos(2*pi*y0(i)/a)-alpha/4*cos(4*pi*y0(i)/a));

	Un_d = 0;

	for j=1:Nsteps
		y01(j) = y0(i)+(j-1)*Dy(i);
		y02(j) = y0(i)+j*Dy(i);

		Gamma_y01(j) = (Gamma_c+Gamma_0)/2+(Gamma_c-Gamma_0)/2*(alpha/4+cos(2*pi*y01(j)/a)-alpha/4*cos(4*pi*y01(j)/a));
		Gamma_y02(j) = (Gamma_c+Gamma_0)/2+(Gamma_c-Gamma_0)/2*(alpha/4+cos(2*pi*y02(j)/a)-alpha/4*cos(4*pi*y02(j)/a));

		A1(j) = sqrt(Gamma_y01(j)^2-(K(i)*tau_P*b*(y01(j)-y0(i)) + Gamma_y0(i))^2);
		A2(j) = sqrt(Gamma_y02(j)^2-(K(i)*tau_P*b*(y02(j)-y0(i)) + Gamma_y0(i))^2);
	Un_d = Un_d + Dy(i)*(A1(j) + A2(j));

	end

	Un = [Un Un_d];

end

%%% Visual check

%lagr1 = [];

%for j = -1:0.05:1

%lagr1 = [lagr1 (Gamma_c-Gamma_0)/2*(cos(2*pi.*j/a) - cos(2*pi*y0(i)/a)) - alpha/8*(Gamma_c - Gamma_0)*(cos(4*pi.*j/a)-cos(4*pi*y0(i)/a))-K(i)*tau_P*b*(j-y0(i))];

%end

%lagr = [lagr; lagr1];

%figure(i); plot([-1:0.05:1]./a,lagr(i,:),'b-',[-2:0.01:2]./a,fun2(-2:0.01:2),'r-',lambda_c(i)/a,0,'r o')

%grid on

