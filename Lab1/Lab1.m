clear

% Дані
data = load("f7.txt");
dataLength = length(data);
deltaF = 0.01;
T = 5;
t = 0:deltaF:T;

% Графік спостережень
figure
plot(t,data), grid

% Дискретне перетворення Фур'є 
fourierFunc = zeros(1,dataLength);
for m = 1:dataLength
    for j = 1:dataLength
        fourierFunc(m) = fourierFunc(m) + 1/dataLength*data(j)*exp(1)^(-1i*2*pi/dataLength*m*j);
    end
end

figure
% delta f і графік перетворення Фур'є для показу екстремумів
deltaF = 1/T;
dataLength = length(t);
plot(abs(fourierFunc)),grid
f = 0:deltaF:round(dataLength/2) * deltaF;

figure
plot(f,abs(fourierFunc(1:round(dataLength/2)+1)))

% Локальні максимуми і частоти
fourierFunc=abs(fourierFunc);
iterator = 0;
extr = zeros(2,1);
for j = 3:round(dataLength/2)-1
if(fourierFunc(j) > fourierFunc(j+1) && fourierFunc(j) > fourierFunc(j-1) && abs(fourierFunc(j)-fourierFunc(j+1)) > 1)
iterator = iterator + 1;
extr(iterator) = j*deltaF;
end
end

% Будуємо та розв'язуємо систему рівнянь, щоб знайти коефіцієнти при частотах
fSin = sin(2*pi*extr(1)*t);
% Розв’язуємо за допомогою матричного методу
A = [sum(t.^6), sum(t.^5), sum(t.^4), sum(fSin.*t.^3), sum(t.^3);
sum(t.^5), sum(t.^4), sum(t.^3), sum(fSin.*t.^2), sum(t.^2);
sum(t.^4), sum(t.^3), sum(t.^2), sum(fSin.*t),    sum(t);
sum(fSin.*t.^3), sum(fSin.*t.^2), sum(fSin.*t), sum(fSin.*fSin), sum(dataLength*fSin);
sum(t.^3), sum(t.^2), sum(t), sum(dataLength*fSin), dataLength];

c = [sum(data.*t.^3), sum(data.*t.^2), sum(data.*t), sum(data.*fSin),  sum(data)];

a = INV(A)*c';
temp = a';

% Отримана апроксимуюча функція
aproxF = a(1).*t.^3 + a(2).*t.^2 + a(3).*t + a(4).*fSin +a(5);

% Графік апроксимуючої функції
figure
plot(t, aproxF), grid

% Середньоквадратична похибка
deviation = immse(data,aproxF);
