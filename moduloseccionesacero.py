import numpy as np
import sympy as sp



#Formulas obtenidas de manual ICHA
E = 200000 #MPa
G = 77200 #MPa

class Perfilhsoldado:
    def __init__(self, h, bf, tf, tw):
        self.h = h
        self.bf = bf
        self.tf = tf
        self.tw = tw
        self.d = self.h + 2*self.tf

    def A(self):
        a = 2 * self.bf * self.tf + self.h * self.tw
        return a

    #Momentos de inercia
    def Ix(self):
        ix = (self.bf * self.d**3 - (self.bf - self.tw)*self.h**3)/12
        return ix   
    
    def Iy(self):
        iy = (2 * self.tf * self.bf**3 + self.h * self.tw**3)/12
        return iy
    
    def rx(self):
        Ix = Perfilhsoldado.Ix(self)
        A = Perfilhsoldado.A(self)
        rx = np.sqrt(Ix / A)
        return rx
    
    def ry(self):
        Iy = Perfilhsoldado.Iy(self)
        A = Perfilhsoldado.A(self)
        ry = np.sqrt(Iy / A)
        return ry
    
    def Sx(self):
        Ix = Perfilhsoldado.Ix(self)
        sx = Ix/(self.d/2)
        return sx
    
    def Sy(self):
        Iy = Perfilhsoldado.Iy(self)
        sy = Iy/(self.bf/2)
        return sy
    
    #Módulos plásticos
    def Zx(self):
        zx = self.bf * self.tf * (self.h + self.tf) + self.tw * self.h**2 / 4
        return zx
    
    def Zy(self):
        zy = self.tf * self.bf**2 /2 + self.h * self.tw**2 /4
        return zy
    
    #Propiedades flexo torsionales
    def J(self):
        j = (2 * self.bf * self.tf**3 + (self.h + self.tf) * self.tw**3) /3
        return j
    
    def Cw(self):
        cw = self.tf * self.bf**3 * (self.h + self.tf)**2 /24
        return cw
    
    def ia(self):
        Iy = Perfilhsoldado.Iy(self)
        Sx = Perfilhsoldado.Sx(self)
        ia = np.sqrt(self.d * Iy / (2* Sx))
        return ia
    
    def it(self):
        it = self.bf * self.tf / self.d
        return it
    
    def X1(self):
        Sx = Perfilhsoldado.Sx(self)
        J = Perfilhsoldado.J(self)
        A = Perfilhsoldado.A(self)
        x1 = np.pi/Sx * np.sqrt(E * G * J * A /2)
        return x1
    
    def X2(self):
        Cw = Perfilhsoldado.Cw(self)
        Iy = Perfilhsoldado.Iy(self)
        Sx = Perfilhsoldado.Sx(self)
        J = Perfilhsoldado.J(self)
        x2 = 4* Cw / Iy * (Sx / G / J)**2
        return x2

class Perfiltsoldado:
    def __init__(self, h, bf, tf, tw):
        self.h = h
        self.bf = bf
        self.tf = tf
        self.tw = tw
        self.d = self.h + self.tf
    
    def A(self):
        a = self.bf * self.tf + self.h * self.tw
        return a
    
    def y(self):
        A = Perfiltsoldado.A(self)
        y = self.tf / 2 + self.h * self.d * self.tw / (2*A)
        return y

    #Centro plástico
    def yp(self):
        area = Perfiltsoldado.A(self)
        if area/2 >= self.bf * self.tf:
            yp = (area/2 - self.bf * self.tf) / self.tw + self.tf
        else:
            yp = area/(2*self.bf)
        return yp
    
    #Momentos de inercia
    def Ix(self):
        y = Perfiltsoldado.y(self)
        ix = self.bf * self.tf**3 / 12 + self.bf * self.tf * (y - self.tf/2)**2 + self.tw * self.h**3 / 12 + self.tw * self.h * (self.h/2 + self.tf - y)**2
        return ix

    def Iy(self):
        iy = (self.tf * self.bf**3 + self.h * self.tw**3)/12
        return iy
    
    def rx(self):
        Ix = Perfiltsoldado.Ix(self)
        A = Perfiltsoldado.A(self)
        rx = np.sqrt(Ix / A)
        return rx
    
    def ry(self):
        Iy = Perfiltsoldado.Iy(self)
        A = Perfiltsoldado.A(self)
        ry = np.sqrt(Iy / A)
        return ry
    
    def Sx_sup(self):
        Ix = Perfiltsoldado.Ix(self)
        y = Perfiltsoldado.y(self)
        sxsup = Ix / y
        return sxsup
    
    def Sx_inf(self):
        Ix = Perfiltsoldado.Ix(self)
        y = Perfiltsoldado.y(self)
        sxinf = Ix / (self.d - y)
        return sxinf
    
    def Sy(self):
        Iy = Perfiltsoldado.Iy(self)
        sy = Iy/(self.bf/2)
        return sy
    
    #Módulo plástico
    def Zx(self):
        yp = Perfiltsoldado.yp(self)
        if yp >= self.tf:
            zx = self.bf * self.tf * (yp - self.tf / 2) + self.tw * (yp - self.tf)**2 / 2 + self.tw * (self.d - yp)**2 /2
        else:
            zx = self.bf * (yp**2 + self.tf**2 /2 - yp * self.tf) + self.h * self.tw * (self.h/2 + self.tf - yp)
        return zx
    
    def Zy(self):
        zy = (self.tf * self.bf**2 + self.h * self.tw**2)/4
        return zy
    
    #Propiedades flexo torsionales
    def J(self):
        j = (self.bf * self.tf**3 + (self.h + self.tf/2)* self.tw**3) / 3
        return j
    
    def Cw(self):
        cw = (self.tf**3 * self.bf**3 / 4 + self.tw**3 * (self.h + self.tf/2)**3)/36
        return cw
    
    def j(self):
        y = Perfiltsoldado.y(self)
        Ix = Perfiltsoldado.Ix(self)
        j = (((self.d - y)**4 - (y - self.tf/2)**4) * self.tw/4 - self.bf * self.tf * (y - self.tf/2)*((y - self.tf/2)**2 + self.bf**2/12)) / (2*Ix) + (y - self.tf/2)
        return j
    
    def r_o(self):
        y = Perfiltsoldado.y(self)
        Ix = Perfiltsoldado.Ix(self)
        Iy = Perfiltsoldado.Iy(self)
        A = Perfiltsoldado.A(self)
        r_o = ((y - self.tf/2)**2 + (Ix + Iy) / A)**0.5
        return r_o
    
    def H(self):
        y = Perfiltsoldado.y(self)
        r_o = Perfiltsoldado.r_o(self)
        h = 1 - ((y - self.tf/2) / r_o)**2
        return h

class Perfilhlaminado:
    def __init__(self, d, bf, tf, tw, r):
        self.d = d
        self.bf = bf
        self.tf = tf
        self.tw = tw
        self.r = r

    def A(self):
        a = 2 * self.bf * self.tf + self.tw * (self.d - 2*self.tf) + (2*self.r)**2 - np.pi * self.r**2
        return a
    
    #Momentos de inercia
    def Ix(self):
        ix = (self.bf * self.d**3 - (self.bf - self.tw) * (self.d - 2*self.tf)**3)/12 + 0.8584*self.r**2 * (self.d/2 - self.tf - 0.2234*self.r)**2 + 0.0302*self.r**4
        return ix
    
    def Iy(self):
        iy = self.tf * self.bf**3 /6 + (self.d - 2*self.tf)*self.tw**3 + (0.8584*self.r**2) * (self.tw/2 + 0.2234*self.r)**2 + 0.0302*self.r**4
        return iy

    def Sx(self):
        sx = Perfilhlaminado.Ix(self) / (self.d/2)
        return sx
    
    def Sy(self):
        sy = Perfilhlaminado.Iy(self) / (self.bf/2)
        return sy
    
    def rx(self):
        Ix = Perfilhlaminado.Ix(self)
        A = Perfilhlaminado.A(self)
        rx = np.sqrt(Ix / A)
        return rx
    
    def ry(self):
        Iy = Perfilhlaminado.Iy(self)
        A = Perfilhlaminado.A(self)
        ry = np.sqrt(Iy / A)
        return ry
    
    #Módulo plástico
    def Zx(self):
        zx = self.bf * self.tf * (self.d - self.tf) + self.tw * (self.d/2 - self.tf)**2 + 0.8584*self.r**2 * (self.d/2 - self.tf - 0.2234*self.r)
        return zx
    
    def Zy(self):
        zy = self.tf * self.bf**2 /2 + (self.d - 2*self.tf) * self.tw**2 /4 + 0.8584*self.r**2 * (self.tw /2 + 0.2234*self.r)
        return zy
    
    #Propiedades flexo torsionales
    def D(self): #Diametro interno encuentro alma-ala
        d = (self.tf**2 + self.tw**2 /4 + 0.2929*self.r * (self.tw + 2*self.tf) + 0.1716*self.r**2) / (self.tf + 0.2929*self.r)
        return d
    
    def alpha(self):
        a = (0.15 + 0.1*self.r / self.tf) * self.tw / self.tf
        return a
    
    def J(self):
        D = Perfilhlaminado.D(self)
        alpha = Perfilhlaminado.alpha(self)
        J = 2 * self.bf * self.tf**3 * (1/3 - 0.21*self.tf * (1 - self.tf**4 / (12 * self.bf**4))/self.bf) + (self.d - 2*self.tf) * self.tw**3 /3 + 2*alpha*D**4
        return J
    
    def Cw(self):
        Iy = Perfilhlaminado.Iy(self)
        cw = Iy * (self.d - self.tf)**2 /4
        return cw

    def ia(self):
        Iy = Perfilhlaminado.Iy(self)
        Sx = Perfilhlaminado.Sx(self)
        ia = np.sqrt(self.d * Iy / (2*Sx))
        return ia

    def it(self):
        it = self.bf * self.tf / self.d
        return it

    def X1(self):
        Sx = Perfilhlaminado.Sx(self)
        J = Perfilhlaminado.J(self)
        A = Perfilhlaminado.A(self)
        X1 = np.pi / Sx * np.sqrt(E * G * J * A /2)
        return X1
    
    def X2(self):
        Sx = Perfilhlaminado.Sx(self)
        Iy = Perfilhlaminado.Iy(self)
        Cw = Perfilhlaminado.Cw(self)
        J = Perfilhlaminado.J(self)
        X2 = 4* Cw / Iy * (Sx / G / J)**2
        return X2

    #104 747

class Perfilllaminado:
    def __init__(self, a, t, R, R1):
        self.a = a
        self.t = t
        self.R = R
        self.R1 = R1
        self.a_aux = a - t - R1
        self.b_aux = t - R1

    def A(self):
        a = self.t * (2 * self.a - self.t) + 0.2146 * (self.R**2 - 2*self.R1**2)
        return a
    
    #Centro de gravedad x e y
    def cg_xy(self):
        A = Perfilllaminado.A(self)
        cg = (6 * self.t * (self.a * (self.a + self.t) - self.t**2) + self.R1**2 * (1.1504*self.R1 - 2.5752 * (self.a + self.t)) + self.R**2 * (2.5752 * self.t + 0.5752 * self.R)) / (12* A)
        return cg
    
    #Centro plástico
    def Xp(self):
        xp = sp.symbols("xp")
        eq = (xp - self.b_aux) * sp.sqrt(self.R1**2 - (xp - self.b_aux)**2) + 2*xp * (self.t + self.a_aux) - self.R1**2 * sp.acos((xp - self.b_aux)/self.R1) - self.t**2 - 0.2146*self.R**2 - 2 * self.t * self.a_aux
        inicial = self.b_aux
        
        if abs(inicial - self.b_aux) >= self.R1:
            raise ValueError("Valor inicial fuera de dominio geométrico al intentar calcular xp e yp")
        
        sol = sp.nsolve(eq, xp, inicial) 
        return sol
    
    #Momentos de inercia
    def I(self):
        y = Perfilllaminado.cg_xy(self)
        inercia = ((self.a_aux * self.t**3 + self.a_aux**3 * self.t + self.t**4 + self.R1 * self.b_aux**3 + self.R1**3 * self.b_aux)/12
        + self.t * (y - self.t/2)**2 * (self.a_aux + self.t) + 0.0075 * self.R**4
        + self.R1 * self.b_aux * ((y - self.b_aux/2)**2 + (self.a_aux + self.t + self.R1/2 - y)**2)
        + 0.7854 * self.R1**2 * ((y - self.b_aux - 0.4244 * self.R1)**2 + (self.a_aux + self.t + 0.4244*self.R1 - y)**2)
        + 0.2146*self.R**2 * (y - self.t - 0.2234 * self.R)**2 + self.a_aux * self.t * (y - self.t - self.a_aux/2)**2)
        return inercia
    
    def Ixy(self):
        x = Perfilllaminado.cg_xy(self)
        ixy = (self.t * (self.t/2 - x) * (self.a**2 - 2*self.a * x + self.t * x - self.t**2/2) 
        - 0.1065*(self.R**4 /24 - self.R1**4 /12) 
        + 0.2146*self.R**2 * (x - self.t - 0.2234*self.R)**2
        - 0.4292*self.R1**2 * (self.a - x - 0.2234*self.R1) * (self.t - x - 0.2234*self.R1))
        return ixy

    def Iu(self):
        return Perfilllaminado.I(self) - Perfilllaminado.Ixy(self)
    
    def ru(self):
        Iu = Perfilllaminado.Iu(self)
        A = Perfilllaminado.A(self)
        ru = np.sqrt(Iu / A)
        return ru
    
    def Iv(self):
        return Perfilllaminado.I(self) + Perfilllaminado.Ixy(self)
    
    def rv(self):
        Iv = Perfilllaminado.Iv(self)
        A = Perfilllaminado.A(self)
        rv = np.sqrt(Iv / A)
        return rv
    
    #Módulo plástico
    def Z(self):
        Xp = Perfilllaminado.Xp(self)
        z = (self.a * (self.t - Xp)**2 + self.t * (self.a**2 - self.t**2 + 2*self.t * Xp - self.a * self.t)/2 
        + self.R**2 * (2.5752 * (self.t - Xp) + 0.5752*self.R)/12 - 0.2146*self.R1**2 * (self.a - self.t))
        return z

    #Propiedades flexo torsionales
    def J(self):
        D = 0.8284*self.t + 0.2426*self.R
        alpha = 0.07 + 0.076*self.R / self.t
        J = (self.a * self.t**3 * (1/3 - 0.21*self.t/self.a * (1- self.t**4 / (12*self.a**4))) 
        + self.t**3 * (self.a - self.t) * (1/3 - 0.105 * self.t / (self.a - self.t) * (1 - self.t**4 / (192 * (self.a - self.t)**4)))) + alpha * D**4
        return J

    def Cw(self):
        cw = self.t**3 * (self.a - self.t/2)**3 /18
        return cw
    
    def x0(self):
        x = Perfilllaminado.cg_xy(self)
        x0 = (x - self.t/2) * np.sqrt(2)
        return x0
    
    def j(self):
        x0 = Perfilllaminado.x0(self)
        Iv = Perfilllaminado.Iv(self)
        j = np.sqrt(2) * self.t * (self.a - self.t/2)**4 / (48 * Iv) + x0
        return j

    def r0(self):
        x0 = Perfilllaminado.x0(self)
        Ix = Perfilllaminado.I(self)
        A = Perfilllaminado.A(self)
        r0 = (x0**2 + 2 * Ix / A)**0.5
        return r0
    
    def H(self):
        x0 = Perfilllaminado.x0(self)
        r0 = Perfilllaminado.r0(self)
        h = 1 - (x0 / r0)**2
        return h

class Perfilcplegado:
    def __init__(self, D, B, t, R):
        self.D = D
        self.B = B
        self.t = t
        self.R = R
        self.r = R + t/2
        self.u = np.pi*self.r/2
        self.a = D - 2*(t + R)
        self._a = D - t
        self.b = B - t - R
        self._b = B - t/2

    def A(self):
        a = self.t * (self.a + 2*(self.b + self.u))
        return a
    
    def cg(self):
        A = Perfilcplegado.A(self)
        x = self.t * (self.a * self.t /2 + self.b**2 + (self.b + self.u)*(2*self.r + self.t) - 2*self.r**2 - self.t**2/6) /A
        return x
    
    #Centro plástico
    def Xp(self):
        A = Perfilcplegado.A(self)
        A1 = (self.r + self.t/2)**2 * float(sp.atan(self.t / self.r)) - self.R * self.r/2
        if self.b * self.t >= A/4:
            xp = self.b/2 + 0.2146*self.r + self.t/2 - self.a/4
            print(f"Caso 1. (Xp = {xp}) >= {self.R + self.t}, eje en el tramo recto del ala")
            return float(xp)
        elif A1 < (A/4 - self.a * self.t/2):
            theta = (A/2 - self.a * self.t) / (2*self.r * self.t)
            xp = self.t/2 + self.r * (1 - float(sp.cos(theta)))
            print(f"{self.t} <= (Xp = {xp}) < {self.R + self.t}, eje en el codo")
            return float(xp)
        else:
            xp = sp.symbols("xp")
            theta2 = sp.atan(sp.sqrt(2*xp * (self.r + self.t/2) - xp**2) / (self.r + self.t/2 - xp))
            eq = sp.nsolve(self.a * xp + (self.r + self.t/2)**2 * (theta2 - 0.5*sp.sin(2*theta2)) - A/2, xp, 1)
            print(f'(Xp = {eq}) < {self.t}, eje en el alma')
            return float(eq)

    #Momentos de inercia
    def Ix(self):
        ix = 2*self.t * (0.0417 * self.a**3 + self.b * (self.a/2 + self.r)**2 + self.u * (self.a/2 + 0.637*self.r)**2 + 0.149 * self.r**3)
        return ix
    
    def Iy(self):
        A = Perfilcplegado.A(self)
        x = Perfilcplegado.cg(self)
        iy = 2*self.t * (0.0833*self.b**3 + self.b * (self.b/2 + self.r)**2 + 0.356*self.r**3) - A*(x - self.t/2)**2
        return iy
    
    def Sx(self):
        Ix = Perfilcplegado.Ix(self)
        return Ix / (self.D/2)

    def Sy_der(self):
        x = Perfilcplegado.cg(self)
        Iy = Perfilcplegado.Iy(self)
        sy_der = Iy / (self.B - x)
        return sy_der

    def Sy_izq(self):
        x = Perfilcplegado.cg(self)
        Iy = Perfilcplegado.Iy(self)
        sy_izq = Iy / (x)
        return sy_izq

    def rx(self):
        A = Perfilcplegado.A(self)
        Ix = Perfilcplegado.Ix(self)
        rx = np.sqrt(Ix / A)
        return rx

    def ry(self):
        A = Perfilcplegado.A(self)
        Iy = Perfilcplegado.Iy(self)
        ry = np.sqrt(Iy / A)
        return ry
    
    #Módulo plástico
    def Zx(self):
        zx = self.t * (self.a**2/4 + self._a * self.b + np.pi * self.r * self.a /2 + 2*self.r**2 + self.t**2 /6)
        return zx
    
    def Zy(self):
        xp = Perfilcplegado.Xp(self)
        A = Perfilcplegado.A(self)
        if xp >= (self.R + self.t):
            zy = (self.t * (self.a * (xp - self.t/2) + float(pi) * self.r * (xp - self.r - self.t/2) 
                + (self.b + self.r + self.t/2 - xp)**2 + (xp - self.r - self.t/2)**2 + 2*self.r**2 + self.t**2 /6))
            print(f"Caso 1: (Xp = {xp} >= {self.R + self.t})")
            return zy
        elif self.t <= xp and xp < (self.R + self.t):
            theta = (A/2 - self.a * self.t) / (2*self.r * self.t)
            zy = (self.t * (xp*(self.a + self.r * (3*theta - float(pi)/2) 
                - 2*self.b) + self.r**2 * (2*float(sp.sin(theta)) 
                - 3*float(pi)*theta /2) + self.r * self.t * (float(pi)/2 - 3*theta)/2 
                + self.b * (2* self.B - self.b) - self.a * self.t/2))
            print(f"Caso 2: {self.t} <= (Xp = {xp} < {self.R + self.t}")
            return zy
        elif xp < self.t:
            theta2 = float(sp.atan(sp.sqrt(2*xp * (self.r + self.t/2) - xp**2) / (self.r + self.t/2 - xp)))
            zy = (self.a * (xp**2 - self.t * xp + self.t**2/2) + 3/8 * xp * (self.r + self.t/2)**2 * (theta2 - 0.5*float(sp.sin(2*theta2))) 
            + 0.5*(float(pi) * self.r * self.t - (self.r + self.t/2)**2 * (theta2 - 0.5*float(sp.sin(2*theta2)))) * (self.r + self.t/2 - xp) + 2*self.b * self.t * (self.B - self.b/2 - xp))
            print(f"Caso 3: (Xp = {xp} < {self.t})")
            return zy

    #Propiedades flexo torsionales
    def m(self):
        m = 3 * self._b**2 / (self._a + 6*self._b)
        return m
    
    def J(self):
        j =self.t**3 * (self.a + 2*self.b + 2*self.u)/3
        return j
    
    def Cw(self):
        cw = (self.t * self._a**2 * self._b**2 /12) * (2* self._a**3 * self._b + 3*self._a**2 * self._b**2) / (6*self._a**2 * self._b + self._a**3)
        return cw
    
    def x0(self): #distancia entre el CG y centro de corte CC
        x = Perfilcplegado.cg(self)
        m = Perfilcplegado.m(self)
        x0 = x + m - self.t/2
        return x0
    
    def j(self):
        x = Perfilcplegado.cg(self)
        x0 = Perfilcplegado.x0(self)
        Iy = Perfilcplegado.Iy(self)
        betaw = -(self.t * self._a**3 * (x - self.t/2)/12 + self.t * self._a * (x - self.t/2)**3)
        betaf = self.t * ((self._b - x + self.t/2)**4 - (x - self.t/2)**4)/2 + self.t * self._a**2 * ((self._b - x + self.t/2)**2 - (x - self.t/2)**2)/4
        j = x0 + (betaw + betaf) / (2* Iy)
        return j
    
    def r0(self):
        x0 = Perfilcplegado.x0(self)
        Ix = Perfilcplegado.Ix(self)
        Iy = Perfilcplegado.Iy(self)
        A = Perfilcplegado.A(self)
        r0 = (x0**2 + (Ix + Iy) / A)**0.5
        return r0
    
    def H(self):
        x0 = Perfilcplegado.x0(self)
        r0 = Perfilcplegado.r0(self)
        H = 1 - (x0/r0)**2
        return H
    
    def ia(self):
        Iy = Perfilcplegado.Iy(self)
        Sx = Perfilcplegado.Sx(self)
        ia = np.sqrt(self.D * Iy / 2 / Sx)
        return ia
    
    def it(self):
        it = self.B * self.t / self.D
        return it
    
    def X1(self):
        J = Perfilcplegado.J(self)
        A = Perfilcplegado.A(self)
        Sx = Perfilcplegado.Sx(self)
        x1 = np.pi / Sx * np.sqrt(E * G * J * A /2)
        return x1

    def X2(self):
        Cw = Perfilcplegado.Cw(self)
        Iy = Perfilcplegado.Iy(self)
        Sx = Perfilcplegado.Sx(self)
        J = Perfilcplegado.J(self)
        x2 = 4* Cw / Iy * (Sx /G /J)**2
        return x2

class Perfilcaplegado:
    def __init__(self, D, B, d, t, R):
        self.D = D
        self.B = B
        self.d = D
        self.t = t
        self.R = R
        self.r = R + t/2
        self.u = np.pi * self.r/2
        self.a = D - 2*(t + R)
        self._a = D - t
        self.b = B - 2*(t + R)
        self._b = B - t/2
        self.c = d - t - R
        self._c = d - t/2

    def A(self):
        a = self.t * (self.a + 2*self.b + 2*self.c + 4*self.u)
        return a
        
    def cg(self):
        A = Perfilcaplegado.A(self)
        x = (self.a * self.t**2 /2 + 2*self.b * self.t * self.B/2 + 2*(self.c * self.t) * (self.B - self.t/2) + 4*self.u * self.t * self.B/2 ) /A
        return x 
    
    def Xp(self):
        A = Perfilcaplegado.A(self)
        A1 = (self.r + self.t/2)**2 * np.arctan(self.t / self.r) - self.R * self.r /2
        if self.t * (self.b + self.c + self.u) >= A/4:
            xp = (self.B + self.d - self.D/2)/2
            print(f"Caso 1: (Xp = {xp}) >= {self.R + self.t}, eje en el tramo recto del ala")
            return xp
        elif A1 < (A/4 - self.a * self.t/2):
            theta = (A/2 - self.a * self.t) / (2* self.r * self.t)
            xp = self.t/2 + self.r * (1 - np.cos(theta))
            print(f"Caso 2: {self.t} <= (Xp = {xp}) < {self.R + self.t}, eje en el codo")
            return xp
        else:
            try:
                xp = sp.symbols("xp")
                A = Perfilcaplegado.A(self)
                theta2 = sp.atan(sp.sqrt(2*xp * (self.r + self.t/2) - xp**2) / (self.r + self.t/2 - xp))
                eq = sp.nsolve(self.a * xp + (self.r + self.t/2)**2 * (theta2 - 0.5*sp.sin(2*theta2) - A/2), xp, 1 )
                print(f"Caso 3: (Xp = {eq}) < {self.t}")
                return float(eq)
            except TypeError:
                return "-"
    
    def Ix(self):
        ix = 2*self.t * (0.0417*self.a**3 + self.b * (self.a/2 + self.r)**2 + 2*self.u*(self.a/2 + 0.637*self.r)**2 + 0.298*self.r**3 + 0.0833*self.c**3 + self.c * (self.a - self.c)**2 /4)
        return ix
    
    def Sx(self):
        Ix = Perfilcaplegado.Ix(self)
        sx = Ix / (self.D/2)
        return sx

    def rx(self):
        Ix = Perfilcaplegado.Ix(self)
        A = Perfilcaplegado.A(self)
        rx = np.sqrt(Ix / A)
        return rx
    
    def Iy(self):
        A = Perfilcaplegado.A(self)
        x = Perfilcaplegado.cg(self)
        iy = 2*self.t * (0.0833*self.b**3 + self.b * (self.b/2 + self.r)**2 + 0.505*self.r**3 + self.c * (self.b + 2*self.r)**2 + self.u*(self.b + 1.637*self.r)**2) - A*(x - self.t/2)**2
        return iy
    
    def Sy_der(self):
        Iy = Perfilcaplegado.Iy(self)
        cg = Perfilcaplegado.cg(self)
        sy = Iy / (self.B - cg)
        return sy
    
    def Sy_izq(self):
        Iy = Perfilcaplegado.Iy(self)
        cg = Perfilcaplegado.cg(self)
        sy = Iy / (cg)
        return sy

    def ry(self):
        Iy = Perfilcaplegado.Iy(self)
        A = Perfilcaplegado.A(self)
        ry = np.sqrt(Iy / A)
        return ry
    
    #Módulo plástico
    def Zx(self):
        zx = self.t * (self.a**2 /4 + self._a * self.b + np.pi*self.r * self.a + 4*self.r**2 + self.t**2/3 + self.c * self.a - self.c**2)
        return zx

    def Zy(self):
        xp = Perfilcaplegado.Xp(self)
        A = Perfilcaplegado.A(self)
        if isinstance(xp, float):
            if xp >= (self.R + self.t):
                zy = self.t * (self.a * (xp - self.t/2) + np.pi * self.r * self.b + (self.b + self.r + self.t/2 - xp)**2 + (xp - self.r - self.t/2)**2 + 4*self.r**2 + self.t**2/3 + 2*self.c*(self.B - xp - self.t/2))
                print(f"Caso 1: (Xp = {xp}) >= {self.R + self.t}")
                return zy
            elif self.t <= xp and xp < (self.R + self.t):
                theta = (A/2 - self.a * self.t) / (2* self.r * self.t)
                zy = self.t * (xp * (self.a + 3*self.r * (theta - np.pi/2) - 2*(self.b + self.c)) + self.r**2 * (2* np.sin(theta) - 3*theta + 2 - np.pi/2) - self.r * self.t * (np.pi/2 + theta)/2 + self.B * (self.b + np.pi*self.r + 2*self.c) + self.t * (self.t/6 - self.a/2 - self.c))
                print(f"Caso 2: {self.t} <= (Xp = {x}) < {self.R + self.t}")
                return zy
            elif xp < self.t:
                theta2 = sp.atan(sp.sqrt(2*xp * (self.r + self.t/2) - xp**2) / (self.r + self.t/2 - xp))
                zy = self.a * (xp**2 - self.t*xp + self.t**2/2) + (7/8*xp - self.r/2 - self.t/4) * (self.r + self.t/2)**2 * (theta2 - 0.5*sp.sin(2*theta2)) + 3/2*float(sp.pi) * self.r * self.t * (self.r + self.t/2 - xp + 2/3*self.b) + self.b * self.t * (self.B + 2*xp) + 2*self.c * self.t * (self.B - self.t/2 - xp) + 2*self.t * self.r**2 + self.t**3/6
                print(f"Caso 3: (Xp = 3) < {self.t}")
                return zy
        return "-"    



#751 - 72
D = 225 #mm
B = 50 #mm
d = 20 #mm
t = 5 #mm
R = 7.5 #mm


perfilcaplegado = Perfilcaplegado(D, B, d, t, R)
print(f"A: {perfilcaplegado.A()} mm2")
print(f"CG: {perfilcaplegado.cg()} mm")
print(f"Xp: {perfilcaplegado.Xp()} mm")
print(f"Ix: {perfilcaplegado.Ix()/10**6} mm4")
print(f"Iy: {perfilcaplegado.Iy()/10**6} mm4")
print(f"Sx: {perfilcaplegado.Sx()/10**3} mm3")
print(f"Sy_der: {perfilcaplegado.Sy_der()/10**3} mm3")
print(f"Sy_izq: {perfilcaplegado.Sy_izq()/10**3} mm3")
print(f"rx: {perfilcaplegado.rx()} mm")
print(f"ry: {perfilcaplegado.ry()} mm")
print(f"Zx: {perfilcaplegado.Zx()/10**3} mm3")
print(f"Zy: {perfilcaplegado.Zy()/10**3} mm3") #ver cómo dividir por 10^3 solo si el valor es un número
# print(f"m: {perfilcaplegado.m()} mm")
# print(f"J: {perfilcaplegado.J()/10**4} mm4")
# print(f"Cw: {perfilcaplegado.Cw()/10**6} mm6")
# print(f"j: {perfilcaplegado.j()} mm")
# print(f"x0: {perfilcaplegado.x0()} mm")
# print(f"r0: {perfilcaplegado.r0()} mm")
# print(f"H: {perfilcaplegado.H()} mm")
# print(f"ia: {perfilcaplegado.ia()} mm")
# print(f"it: {perfilcaplegado.it()} mm")
# print(f"X1: {perfilcaplegado.X1()} MPa")
# print(f"X2: {perfilcaplegado.X2()*10**8} (1/MPa)^2")


