import modulo_secciones_acero as secciones
from materiales import E, fy
E = E*10.19 #kg/cm2

def clasificar(razon, lim_inf, lim_sup):
    if razon <= lim_inf:
        return "Compacto - No Esbelto"
    elif lim_inf < razon < lim_sup:
        return "No Compacto - No Esbelto"
    elif lim_sup <= razon:
        return "No Compacto - Esbelto"


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
    def __init__(self, b=None, t=None, bfs=None, bfi=None, tfs=None, tfi=None, h=None, d=None, tw=None, D=None):
        if b is not None:
            self.b = b / 2
        elif bfs is not None:
            self.b = bfs/2 #Considerar bfs (ancho de ala superior) como el ancho del ala en compresión        
        self.bfs = bfs
        self.bfi = bfi
        self.t = t
        self.tfs = tfs
        self.tfi = tfi
        self.h = h #altura alma
        self.d = d
        self.tw = tw
        self.D = D

    def caso10(self):
        if self.bfs is None or self.tfs is None:
            raise ValueError("bfs y tfs son obligatorios para el caso 10")
        
        razon = self.b / self.t
        lim_inf = 0.38 * (E / fy)**0.5
        lim_sup = 1 * (E / fy)**0.5
        return clasificar(razon, lim_inf, lim_sup)
        
    def caso11(self):
        razon = self.b / self.tfs
        #Cálculo de kc
        kc = 4 / (self.h / self.tw)**0.5
        kc = 0.35 if kc < 0.35 else kc
        kc = 0.76 if kc > 0.76 else kc
        print(f'kc = {kc}')

        #Cálculo de FL
        Sxt = secciones.Perfilhsoldado(self.h, self.bfs, self.bfi, self.tfs, self.tfi, self.tw).Sx_inf()
        Sxc = secciones.Perfilhsoldado(self.h, self.bfs, self.bfi, self.tfs, self.tfi, self.tw).Sx_sup()
        if Sxt / Sxc >= 0.7:
            FL = 0.7 * fy
        if Sxt / Sxc < 0.7:
            FL = fy * Sxt / Sxc
            FL = 0.5*fy if FL < 0.5*fy else FL
        print(f'Sxt/Sxc = {Sxt/Sxc} => FL = {FL}')

        #Clasificación del ala
        lim_inf = 0.38 * (E / fy)**0.5
        lim_sup = 0.95 * (kc * E / FL)**0.5
        return clasificar(razon, lim_inf, lim_sup)
    
    def caso12(self):
        if self.b is None or self.t is None:
            raise ValueError("b y t son obligatorios para el caso 12")
        
        razon = self.b / self.t
        lim_inf = 0.54 * (E / fy)**0.5
        lim_sup = 0.91 * (E / fy)**0.5
        return clasificar(razon, lim_inf, lim_sup)

    def caso13(self):
        if self.b is None or self.t is None:
            raise ValueError("b y t son obligatorios para el caso 13")
        
        razon = self.b / self.t
        lim_inf = 0.38 * (E / fy)**0.5
        lim_sup = 1 * (E / fy)**0.5
        return clasificar(razon, lim_inf, lim_sup)
    
    def caso14(self):
        if self.d is None or self.t is None:
            raise ValueError("d y t son obligatorios para el caso 14")
        
        razon = self.d / self.t
        lim_inf = 0.84 * (E / fy)**0.5
        lim_sup = 1.03 * (E / fy)**0.5
        return clasificar(razon, lim_inf, lim_sup)
    
    ### ELEMENTOS ATIESADOS ###
    
    def caso15(self):
        if self.h is None or self.tw is None:
            raise ValueError("h y tw son obligatorios para el caso 15")
        
        razon = self.h / self.tw
        lim_inf = 3.76 * (E / fy)**0.5
        lim_sup = 5.7 * (E / fy)**0.5
        return clasificar(razon, lim_inf, lim_sup)
    
    def caso16(self):
        perfil = secciones.Perfilhsoldado(self.h, self.bfs, self.bfi, self.tfs, self.tfi, self.tw)
        Sx_sup = perfil.Sx_sup()
        Zx = perfil.Zx()
        hp = perfil.hp()
        hc = perfil.hc()
        My = fy * Sx_sup
        Mp = fy * Zx

        razon = hc / self.tw
        lim_inf1 = hc / hp * (E / fy)**0.5
        lim_inf2 = (0.54 * Mp / My - 0.09)**2
        lim_inf = lim_inf1 / lim_inf2
        lim_sup = 5.7 * (E / fy)**0.5
        lim_inf = lim_sup if lim_inf > lim_sup else lim_inf
        print(razon)
        print(lim_inf)
        print(lim_sup)
        return clasificar(razon, lim_inf, lim_sup)

    def caso17(self):
        if self.b is None or self.t is None:
            raise ValueError("b y t son obligatorios para el caso 17")
        
        razon = self.b / self.t
        lim_inf = 1.12 * (E / fy)**0.5
        lim_sup = 1.4 * (E / fy)**0.5
        return clasificar(razon, lim_inf, lim_sup)
    
    def caso18(self):
        if self.b is None or self.t is None:
            raise ValueError("b y t son obligatorios para el caso 18")
        
        razon = self.b / self.t
        lim_inf = 1.12 * (E / fy)**0.5
        lim_sup = 1.4 * (E / fy)**0.5
        return clasificar(razon, lim_inf, lim_sup)
    
    def caso19(self):
        if self.h is None or self.t is None:
            raise ValueError("h y t son obligatorios para el caso 19")
        
        razon = self.h / self.t
        lim_inf = 2.42 * (E / fy)**0.5
        lim_sup = 5.7 * (E / fy)**0.5
        return clasificar(razon, lim_inf, lim_sup)
    
    def caso20(self):
        if self.D is None or self.t is None:
            raise ValueError("D y t son obligatorios para el caso 20")
        
        razon = self.D / self.t
        lim_inf = 0.07 * E / fy
        lim_sup = 0.31 * E / fy
        return clasificar(razon, lim_inf, lim_sup)
        

