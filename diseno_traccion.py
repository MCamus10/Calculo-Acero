from factores_de_resistencia import PHI
from materiales import E, fy, fu

def flexion(seccion, metodo="LRFD"):
    Ag = seccion.A() /10**2 #cm2
    # Ae = seccion.A_efectiva()
    phi_f = PHI[metodo]["traccion-fluencia"]
    phi_r = PHI[metodo]["traccion-rotura"]

    #Fluencia en tracción de sección bruta
    Pn = fy * Ag
    phiPn = Pn * phi_f
    return phiPn

def flexion(seccion, metodo="ASD"):
    Ag = seccion.A() /10**2 #cm2
    # Ae = seccion.A_efectiva()
    phi_f = PHI[metodo]["traccion-fluencia"]
    phi_r = PHI[metodo]["traccion-rotura"]

    #Fluencia en tracción de sección bruta
    Pn = fy * Ag
    phiPn = Pn / phi_f
    return phiPn