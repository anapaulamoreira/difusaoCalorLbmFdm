%codigo para solucao  da equacao da difusao de calor em uma placa
%condicoes de contorno Dirichlet

clc
close all
clear all

%% M�todo Lattice Boltzmann

n = 100;    % n�mero de n�s de lattice
m = n + 1;  % n�mero de n�s de lattice + 1
dt=1;       % intervalo de tempo
dx=1;       % comprimento de ligacao

%Inicializando os vetores
f1 = zeros(1,m);
f2 = zeros(1,m);
T = zeros(1,m);
feq = zeros(1,m);
x = zeros(1,m);

x(1) = 0;   % comprimento inicial da placa

% comprimento total da placa
for i=2:m
    x(i) = x(i-1) + dx;
end

csq = dx*dx/(dt*dt);
alpha = 0.25;       % coeficiente de difus�o t�rmica
omega = 1/(alpha/(dt*csq) + 0.5);  % frequ�ncia de colis�o

mstep=200;  % n�mero total de etapas de tempo
twall=1;    % Temperatura em que a superficie esquerda � submetida

%% Condi��o inicial

for i=1:m
    T(i)=0; % Valor inicial da temperatura do dom�nio
    f1(i)=0.5*T(i);
    f2(i)=0.5*T(i);
end

%%  main loop
for kk=1:mstep
% processo de colisao:
    for i=1:m
        T(i)=f1(i)+f2(i);
        feq(i)=0.5*T(i);
% feq1=feq2=feq
        f1(i)=(1 - omega)*f1(i)+ omega*feq(i);
        f2(i)=(1 - omega)*f2(i)+ omega*feq(i);
    end
% processo adveccao:
    for i=2:m-2
    f1(m-i) = f1(m-i-1); % f1
    f2(i-1) = f2(i);     % f2
    end
    
% Condi��o de contorno

    f1(1)=twall-f2(1);   % temperatura constante, x = 0   
    f1(m)=f1(m-1);       % adiab�tica, x=L
    f2(m)=f2(m-1);       % adiab�tica, x=L
end


%% Diferen�as Finitas

fo = zeros(1,m);
f = zeros(1,m);
dxd = 1.0;
dtd = 0.500;    % passo de tempo
mstepd = 400;   % numero total de iteracoes

fo(1) = 1.0;    % condicao inicial para o valor antigo de f em x = 0.
f(1) = 1.0;     % condicao inicial para o novo valor de f em x = 0.
fo(m) = fo(m-1);% condicao inicial para o valor antigo de f em x=L
f(m) = f(m-1);  % condicao inicial para o novo valor de f em x=L

for kk=1:mstepd
    for i=2:m-1
        f(i)=fo(i)+dtd*alpha*(fo(i+1)-2.*fo(i)+fo(i-1))/(dxd*dxd);
    end
    
    for i=2:m-1
        fo(i)=f(i); % atualizando o valor antigo
    end 
    fo(m)=f(m-1);   % atualizando a condi��o de contorno em x = L
    end 

% plot dos resultados em um grafico
plot (x, T, 'o', x, f, '+', 'MarkerSize', 10, 'LineWidth', 1)
pause(0.1)
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);
set(gca,'FontSize',15)
xlabel('x');
ylabel('T');
legend('LBM','FDM')
ylim([-0.1 1.1])
print('DirichletLBMFDM','-dpng')