from materiales import E, fy
E = E*10.19 #kg/cm2

class Elemento_en_compresion_por_compresion:
    def __init__(self, b=None, t=None, tw=None, h=None, d=None, D=None):
        self.b = b #Ancho ala
        self.t = t #Espesor ala
        self.h = h #Alto alma
        self.tw = tw #Espesor alma
        self.d = d #Alto sección T
        self.D = D #Diámetro sección tubular
        self.E = 2100000 #kg/cm2
        self.fy = 4200 #kg/cm2

    def caso1(self):
        if self.b is None or self.t is None:
            raise ValueError("Para Caso 1 se requieren parámetros 'b' y 't'")
        razon = self.b / self.t
        limite = 0.56 * (self.E / self.fy)**0.5
        if razon >= limite:
            print("Elemento Ala esbelto")
            return True #true = esbelto
        elif razon < limite:
            print("Elemento Ala no esbelto")
            return False
        
    def caso2(self):
        if self.tw is None or self.h is None or self.b or self.t:
            raise ValueError("Para Caso 2 se requieren los parámetros 'd', 't', 'tw' y 'h'")
        razon = self.b / self.t
        kc = 4 / (self.h / self.tw)
        kc = max(0.35, min(kc, 0.76))

        limite = 0.64 * (kc * self.E / self.fy)**0.5
        if razon >= limite:
            print("Elemento Ala esbelto")
            return True #true = esbelto
        elif razon < limite:
            print("Elemento Ala no esbelto")
            return False
        
    def caso3(self):
        if self.b is None or self.t is None:
            raise ValueError("Para Caso 3 se requieren parámetros 'b' y 't'")
        razon = self.b / self.t
        limite = 0.45 * (self.E / self.fy)**0.5
        if razon >= limite:
            print("Elemento Ala esbelto")
            return True #true = esbelto
        elif razon < limite:
            print("Elemento Ala no esbelto")
            return False
        
    def caso4(self):
        if self.d is None or self.t is None:
            raise ValueError("Para Caso 4 se requiere el parámetro 'd' y 't'")
        razon = self.d / self.t
        limite = 0.75 * (self.E / self.fy)**0.5
        if razon >= limite:
            print("Elemento Alma T esbelto")
            return True #true = esbelto
        elif razon < limite:
            print("Elemento Alma T no esbelto")
            return False
        
    def caso5(self):
        if self.h is None or self.tw is None:
            raise ValueError("Para Caso 5 se requieren parámetros 'h' y 'tw'")
        razon = self.h / self.tw
        limite = 1.49 * (self.E / self.fy)**0.5
        if razon >= limite:
            print("Elemento Alma esbelto")
            return True #true = esbelto
        elif razon < limite:
            print("Elemento Alma no esbelto")
            return False
        
    def caso6(self):
        if self.b is None or self.t is None:
            raise ValueError("Para Caso 6 se requieren parámetros 'b' y 't'")
        razon = self.b / self.t
        limite = 1.4 * (self.E / self.fy)**0.5
        if razon >= limite:
            print("Elemento Ala esbelto")
            return True #true = esbelto
        elif razon < limite:
            print("Elemento Ala no esbelto")
            return False
        
    def caso7(self):
        if self.b is None or self.t is None:
            raise ValueError("Para Caso 7 se requieren parámetros 'b' y 't'")
        razon = self.b / self.t
        limite = 1.4 * (self.E / self.fy)**0.5
        if razon >= limite:
            print("Elemento Ala esbelto")
            return True #true = esbelto
        elif razon < limite:
            print("Elemento Ala no esbelto")
            return False
            
    def caso8(self):
        if self.b is None or self.t is None:
            raise ValueError("Para Caso 8 se requieren parámetros 'b' y 't'")
        razon = self.b / self.t
        limite = 1.49 * (self.E / self.fy)**0.5
        if razon >= limite:
            print("Elemento Ala esbelto")
            return True #true = esbelto
        elif razon < limite:
            print("Elemento Ala no esbelto")
            return False
        
    def caso9(self):
        if self.D is None or self.t is None:
            raise ValueError("Para Caso 9 se requieren parámetros 'D' y 't'")
        razon = self.D / self.t
        limite = 0.11 * self.E / self.fy
        if razon >= limite:
            print("Elemento Circular esbelto")
            return True #true = esbelto
        elif razon < limite:
            print("Elemento Circular no esbelto")
            return False

class Elemento_en_compresion_por_flexion:
    def __init__(self, b=None, t=None, h=None, d=None, tw=None, D=None, My=None, Mp=None):
        self.b = b
        self.t = t
        self.h = h
        self.d = D
        self.tw = tw
        self.D = D
        self.My = My
        self.Mp = Mp

