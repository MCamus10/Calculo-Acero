import numpy as np

class ColumnaHormigon:
    def __init__(self, b, h, fc, fy, recubrimiento):
        """
        b, h: Dimensiones de la sección (mm)
        fc: Resistencia a compresión del hormigón (MPa)
        fy: Tensión de fluencia del acero (MPa)
        """
        self.b = b
        self.h = h
        self.fc = fc
        self.fy = fy
        self.Es = 200000  # Módulo elástico acero (MPa)
        self.beta1 = self._calcular_beta1()
        
        # Lista para almacenar barras de acero: tuplas (Area, distancia_d)
        self.acero = [] 

    def _calcular_beta1(self):
        """Calcula beta1 según ACI 318"""
        if self.fc <= 28:
            return 0.85
        elif self.fc >= 55:
            return 0.65
        else:
            return 0.85 - 0.05 * (self.fc - 28) / 7

    def agregar_barra(self, area, d):
        """
        area: Área de la barra o capa de barras (mm2)
        d: Distancia desde la fibra superior comprimida al centroide de la barra (mm)
        """
        self.acero.append({'As': area, 'd': d})

    def calcular_fuerza_acero(self, deformacion):
        """
        Modelo elasto-plástico perfecto.
        Retorna tensión (MPa)
        """
        # Ley de Hooke: stress = E * strain
        stress = deformacion * self.Es
        
        # Limitar a fy (fluencia)
        stress = np.clip(stress, -self.fy, self.fy)
        return stress

    def calcular_diagrama(self, num_puntos=100):
        puntos_pn = []
        puntos_mn = []
        
        # Centroide plástico (asumido en h/2 para sección rectangular simétrica)
        centroide = self.h / 2

        # 1. Definir rango de eje neutro 'c'
        # Desde compresión pura (c -> infinito) hasta tensión pura (c -> 0 o negativo)
        # Usamos un rango práctico: desde un c muy grande hasta un c muy pequeño
        c_vals = np.linspace(self.h * 3, -self.h * 0.1, num_puntos)
        
        # Agregamos puntos clave manualmente para precisión
        c_vals = np.append(c_vals, [self.h, 0.001]) # c=h (inicio tracción), c~0 (tensión pura)
        c_vals = np.sort(c_vals)[::-1] # Ordenar de mayor a menor (Compresión -> Tensión)

        for c in c_vals:
            Pn = 0
            Mn = 0
            
            # --- HORMIGÓN ---
            # Si c es muy grande, el bloque a se limita a h
            if c > 0:
                a = self.beta1 * c
                a = min(a, self.h) # El bloque no puede salir de la sección
                
                # Fuerza compresión hormigón (0.85 * f'c * a * b)
                Cc = 0.85 * self.fc * a * self.b
                
                # Momento del hormigón respecto al centroide
                # Brazo = (h/2) - (a/2)
                Mc = Cc * (centroide - a/2)
            else:
                Cc = 0
                Mc = 0

            Pn += Cc
            Mn += Mc

            # --- ACERO ---
            # Hipótesis: Deformación máxima del hormigón ecu = 0.003
            # Por semejanza de triángulos: e_s / (c - d) = ecu / c
            # e_s = 0.003 * (c - d) / c
            
            ecu = 0.003
            
            for barra in self.acero:
                d = barra['d']
                As = barra['As']
                
                # Calcular deformación unitaria en la barra
                # Cuidado con c=0
                if abs(c) < 1e-6: 
                    # Tensión pura uniforme (teórica)
                    es = 0.01 # Un valor alto arbitrario para garantizar fluencia en tracción
                else:
                    es = ecu * (c - d) / c
                
                # Calcular tensión fs
                fs = self.calcular_fuerza_acero(es)
                
                # Fuerza en la barra (Positivo = Compresión, igual que el hormigón)
                Fs = As * fs
                
                # Momento de la barra respecto al centroide
                # Brazo positivo hacia arriba desde el centroide: (h/2 - d)
                # Si la barra está abajo (d > h/2), el brazo es negativo, lo cual es correcto.
                Ms = Fs * (centroide - d)
                
                Pn += Fs
                Mn += Ms
            
            # Guardar punto (convertir a kN y kN-m para legibilidad)
            puntos_pn.append(Pn / 1000)      # kN
            puntos_mn.append(Mn / 1000000)   # kN-m

        return puntos_mn, puntos_pn

# --- EJEMPLO DE USO ---

# 1. Crear columna 30x50 cm, H30 (30 MPa), A630 (420 MPa)
col = ColumnaHormigon(b=300, h=500, fc=30, fy=420, recubrimiento=40)

# 2. Definir armadura (Ej: 3fi22 arriba y 3fi22 abajo)
# Área fi22 aprox 380 mm2
area_barra = 380 
# Capa superior (compresión inicial)
col.agregar_barra(area=area_barra*3, d=50) # recubrimiento + estribo + radio
# Capa inferior (tracción inicial)
col.agregar_barra(area=area_barra*3, d=450) # h - 50

# 3. Calcular puntos
momentos, axiales = col.calcular_diagrama(num_puntos=150)

# 4. Imprimir algunos puntos clave para verificar
print(f"{'Moment (kN-m)':<15} | {'Axial (kN)':<15}")
print("-" * 33)
# Imprimimos puntos seleccionados del array (inicio, medio, fin)
indices = [0, len(axiales)//4, len(axiales)//2, 3*len(axiales)//4, -1]
for i in indices:
    print(f"{momentos[i]:<15.2f} | {axiales[i]:<15.2f}")