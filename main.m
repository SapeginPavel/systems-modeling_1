clear all; % очистка рабочей области
nf = 2; % количество факторов
minf = [1 2]; % минимальные значения факторов
maxf = [2 5]; % максимальные значения факторов

% формирование дробного двухуровневого плана эксперимента
fracplan = fracfact('a k ak'); % передача факторов в функцию
fracplan

N = 2 ^ nf; % число экспериментов

% добавление фиктивного фактора
fictfact = ones(N,1);
X = [fictfact fracplan]';

% формируем исходные данные для проведения эксперимента:
fraceks = zeros(N,nf);
for i=1:nf
    for j=1:N
        % Пересчет плана в реальные значения факторов:
        fraceks(j,i) = minf(i) + (fracplan(j,i) + 1) * (maxf(i) - minf(i)) / 2;
    end
end
fraceks

% тактическое планирование эксперимента
d_m = 0.07; % доверительный интервал
alpha = 0.03; % уровень значимости

% определение t-критического
tkr_alpha = norminv(1 - alpha / 2);

% массив для хранения результатов экспериментов
Y = zeros(1, N);

% цикл по совокупности экспериментов стратегического плана
for j = 1:N
    a = fraceks(j,1);
    k = fraceks(j,2);
    
    % расчет количества повторений эксперимента
    D_tilda = (k^2 * pi^2) / 3;
    NE = round(tkr_alpha^2 * D_tilda / d_m^2);
    
    % цикл статистических испытаний
    u = zeros(1, NE);
    for m = 1:NE
        u(m) = systemeqv(a, k); % имитация функционирования системы
    end
    
    % оценка математического ожидания реакции системы
    Y(j) = mean(u);

    figure;
    hist(u,12);
end

% определение коэффициентов регрессии
C = X * X';
b_ = inv(C) * X * Y';

% формирование зависимости реакции системы
A = minf(1):0.1:maxf(1);
B = minf(2):0.1:maxf(2);
[k, N1] = size(A);
[k, N2] = size(B);

Yc = zeros(N2, N1);
Yo = zeros(N2, N1);

for i = 1:N1
    for j = 1:N2
        an = 2 * (A(i) - minf(1)) / (maxf(1) - minf(1)) - 1;
        bn = 2 * (B(j) - minf(2)) / (maxf(2) - minf(2)) - 1;
        Yc(j, i) = b_(1) + an * b_(2) + bn * b_(3) + an * bn * b_(4);
        Yo(j, i) = A(i); % так как для логистического распределения m = a
    end
end

% отображение зависимостей в 3D графике
[x, y] = meshgrid(A, B);
figure;
subplot(1,2,1), plot3(x, y, Yc),
xlabel('Factor a'), ylabel('Factor k'), zlabel('Yc'),
title('System Output - Experimental'), grid on;
subplot(1,2,2), plot3(x, y, Yo),
xlabel('Factor a'), ylabel('Factor k'), zlabel('Yo'),
title('System Output - Theoretical'), grid on;
