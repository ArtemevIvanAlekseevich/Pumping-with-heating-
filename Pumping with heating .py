import matplotlib.pyplot as plt

po =925.19 ##float(input("введите плотность в кг/м^3\n"))

C =1996.53 ##float(input("введите теплоемкость Дж/кг*К\n"))

e =10000 ##float(input("введите максимальный расход в м^3/ч\n"))
e /= 3600

b = 0 / 3600

viz =562.18 ##float(input("введите вязкость при известной Т в мм^2/с\n"))
viz = viz / (10 ** 6)

u =0.1039 ##float(input("введите коэфициент крутизны вязкограммы в 1/К\n"))

Tiz =313 ##float(input("введите Т для вязкости известной в К\n"))

Tн =327.6 ##float(input("введите Тн в К\n"))

T0 =253 ##float(input("введите Т0 в К\n"))

d =315 ##float(input("введите внутренний диаметр в мм\n"))
d /=1000

k =0.1 ##float(input("введите шероховатость в мм\n"))
k /=1000

K1 =0.79 ##float(input("введите полный коэфициент теплоотдачи в ,,,\n"))

x = []
with open("C:\\Users\\user\\Desktop\\Programs\\python\\all\\x.txt") as f:
    for line in f:
        x.append([float(x) for x in line.split()])

z = []
with open("C:\\Users\\user\\Desktop\\Programs\\python\\all\\z.txt") as f:
    for line in f:
        z.append([float(x) for x in line.split()])

n = int (len(x))
Lt = 0
H = []
T = []


for i in range (n):
    H.append (0.0)
    T.append (0.0)
v = 0

T [0] = Tн

Q = 0
t = 0
H0 = 2 * (8.055556 - 0.02778 * 3600 * Q) / 919.415 / 9.81 * 10 ** 6

H [ n - 1 ] = 30 + 109

while abs ((H0 + 151 - H [0]) / (H0 + 151)) * 100 > 0.1:
    Lt = 0
    Q = (b + e) / 2
   

    H0 = 2 * (8.055556 - 0.02778 * 3600 * Q) / 919.415 / 9.81 * 10 ** 6
    
    w = 4 * Q / (3.14 * d ** 2)
    
    for i in range(n-1):
        
        v = viz * 2.7182818284 ** (-u * (T [i] - Tiz))
        
        Re = w *d / ( v)
        if Re < 2000:
            λ = 64 / Re
        elif Re < 10000:
            λ = 64 / Re * (2.7182818284 ** (4.64 - 0.002 * Re)) + 0.3164 / (Re ** 0.25) * (1 - 2.7182818284 ** (4.64 - 0.002 * Re))
        elif Re < 10 * d / k:
            λ = 0.3164 / (Re ** 0.25)
        elif Re < 500 * d / k:
            λ = 0.11 * (68 / Re + k / d) ** 0.25
        else:
            λ = 0.11 * (k/d) ** 0.25
        T [i + 1] = T [i] - (4 * K1 * (T [i] - T0) / (C * po * w * d) - 1 * λ * w * w /(2 * d * C)) * 100
    
    for j in range(n - 1):
        
        v = viz * 2.7182818284 ** (-u * (T [n - j - 1] - Tiz))
        Re = 4 * Q / (3.14 * d * v)
        if Re < 2000:
            λ = 64 / Re
        elif Re < 10 * d / k:
            λ = 0.3164 / (Re ** 0.25)
            Lt += 100
        elif Re < 500 * d / k:
            λ = 0.11 * (68 / Re + k / d) ** 0.25
            Lt += 100
        else:
            λ = 0.11 * (k/d) ** 0.25
            Lt += 100
        
        H [n - j - 2] = H [n - j - 1] + λ * (w ** 2) / (d * 2 * 9.81) * 100 
        
    if H [0] > H0 + 151:
        e = Q
    else:
        b = Q

f = open("C:\\Users\\user\\Desktop\\Programs\\python\\all\\T.txt", "w")
for i in range (n):
    f.write(str (round (T[i],2)))
    f.write('\n')
f.close()

f = open("C:\\Users\\user\\Desktop\\Programs\\python\\all\\H.txt", "w")
for i in range (n):
    f.write(str (round (H[i],2)))
    f.write('\n')
f.close()


print ('Расход ',round(Q * 3600, 2)  , 'м^3/ч')
print ('Напор ',round (H[0]-151, 2), 'м')
print ('Температура конечная',round (T[n-1], 2), 'К')
print ('Длина турбулентного ',Lt , ' м')
print ('Длина Ламинарного ',(n - 1) * 100 - Lt , ' м')

plt.plot(x, z, label = "Профиль трассы")
plt.plot(x, H, label = "Линия гидроуклона")
plt.plot(x, T, label = "температура")

# naming the x axis
plt.xlabel('x - axis')
# naming the y axis
plt.ylabel('y - axis')
# giving a title to my graph
plt.title('Линия гидроуклона')
  
# show a legend on the plot
plt.legend()
  
# function to show the plot
plt.show()

