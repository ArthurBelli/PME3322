# Arthur Pedroso Porto Belli 11804608
# Lucas Kairuz Martins 11805377
# Maria Eduarda Dall Orto de Araujo 11820301 


#constantes do problema
D = 0.05 #diâmetro m
RUG = 0.001 #rugosidade relativa m
LA = 50 #comprimento do trecho A m
LB = 100 #comprimento do trecho B m
LC = 200 #comprimento do trecho C m
LEQ = 50 #comprimento equivalente da válvula C LEQ = K_s*D/f
RHO_H2O = 1000 #massa específica H2O kg/m^3 
NU_H2O = 10e-6 #viscosidade cinemática H2O m^2/s
Z1 = (8+8+8)/3 #cota z1 m
Z2 = Z1 + (0+0+2) + 2 #cota z2 m
Z3 = Z2 + (4+5+0) + 2 #cota z3 m
HM = Z3 + (6+3+3) + 10 #cota manométrica m
G = 10 #aceleração gravitacional m/s^2