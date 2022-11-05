# Arthur Pedroso Porto Belli 11804608
# Lucas Kairuz Martins 11805377
# Maria Eduarda Dall Orto de Araujo 11820301 

from math import *

#constantes do problema
D = 0.05 #diâmetro m
RUG = 0.001 #rugosidade relativa (epsilon/D)
LA = 50 #comprimento do trecho A m
LB = 100 #comprimento do trecho B m
LC = 200 #comprimento do trecho C m
LEQ = 50 #comprimento equivalente da válvula C LEQ = K_s*D/f
RHO_H2O = 1e3 #massa específica H2O kg/m^3 
NU_H2O = 1e-6 #viscosidade cinemática H2O m^2/s
Z1 = (8+8+8)/3 #cota de A m 
Z2 = Z1 + (0+0+2) + 2 #cota de B m 
Z3 = Z2 + (4+5+0) + 2 #cota de C m 
HM = Z3 + (6+3+3) + 10 #cota manométrica m 
G = 10 #aceleração gravitacional m/s^2
REY = D/NU_H2O

'''
Deseja-se saber:
- a carga Hj na junção;
- os coefs de perda de carga distribuída fa, fb, fc ;
- as vazões Qa, Qb, Qc.
'''

#método recursivo para obtenção dos coefs de perda de carga distribuída
def Colebrook(f):
    return (-2*log10(RUG/3.7 + 2.51/((REY)*sqrt(f))))**(-2)
    

def Haaland(reynolds):
    return (-1.8*log10(6.9/reynolds + (RUG/3.7)**(1.11)))**(-2)

def solve_Hj(Hj_ini):
    fcVc2 = (Hj_ini-Z3)*2*G*D/(LC+LEQ) #eq (1)
    fc = Colebrook(fcVc2) #eq (3)
    Vc = sqrt(fcVc2/fc)
    Qc = Vc*pi*D**2/4

    fbVb2 = (Hj_ini-Z2)*2*G*D/(LB) #eq (2) 
    fb = Colebrook(fbVb2) #eq (3)
    Vb = sqrt(fbVb2/fb) #eq (4)
    Qb = Vb*pi*D**2/4

    Va = Vb + Vc #eq (5)
    Qa = Va*pi*D**2/4
    fa = Haaland(Va*REY) #eq (6)

    Hj = Z1 + HM - fa*(LA*(Va**2))/(D*2*G) #eq (7)
    HJ = (Hj + Hj_ini)/2

    if abs(HJ-Hj_ini) < 0.0001: #critério de convergência adotado
        return Hj, Qa, Qb, Qc, fa, fb, fc
    else:
        return solve_Hj(HJ)

def main():
    initial_hj = float(input("Chute inicial de Hj: "))
    Hj, Qa, Qb, Qc, fa, fb, fc = solve_Hj(initial_hj)
    print(f"Hj = {Hj}")
    print(f"fa = {fa}")
    print(f"fb = {fb}")
    print(f"fc = {fc}")
    print(f"Qa = {Qa}")
    print(f"Qb = {Qb}")
    print(f"Qc = {Qc}")

main()