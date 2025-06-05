$offsymxref offsymlist
option limcol =0, limrow=0;
option solprint=off;
option decimals=8;

Sets 
  t         "Tidsperioder"                                  /2024*2054/
  t0(t)     "Første tidsperiode"
  t1(t)     "Første endogene tidsperiode (Basisår)" 
  tT(t)     "Sidste tidsperiode T"
  txT(t)    "Undtagen sidste tidsperiode T"
  
  a         "Bygningsaldre"                                 /0*5/
  a0(a)     "Alder for nye bygninger"                       
  aA(a)     "Højeste bygningsalder A"                         
  axA(a)    "Undtagen højeste bygningsalder A"
;  
 
t0(t)     = yes$(ord(t) = 1);
t1(t)     = yes$(ord(t) = 2);
tT(t)     = yes$(ord(t) = card(t));
txT(t)    = yes$(ord(t) < card(t));

a0(a)     = yes$(ord(a) = 1);
aA(a)     = yes$(ord(a) = card(a));
axA(a)    = yes$(ord(a) < card(a));


Variables
* Til 7-ligningmodellen (se separat fil)
    K(a,t)    "Bygninger"
    K_new(t)  "Nye bygninger"
    K_old(t)   "Gamle bygninger"
    
    p_uc(a,t) "User cost for bygninger"
    p(a,t)    "Markedspris for bygninger"
    p_old(t)  "Markedspris for gamle bygninger"
    p_new(t)  "Markedspris for nye bygninger"
    p_U(t)    "CES-prisindeks for nye og gamle bygninger"

    mu(a)     "Vægt"
    gamma_new "Vægt"
    gamma_old "Vægt"
    

*  Til 15-ligningsmodellen
    W_hat(t)  "Samlet formue"
*    W_c(t)      "Finansiel formue"

    u(t)      "Nytte"
    C_big(t)  "Byggeomkostninger"
    Y_H(t)      "Løbende indtægter"
    G_bar(t)  "Nye arealer til nye bygninger"
    N_old(t)  "Nedrivning af gamle bygninger"
    
;


variables
    E         "Substitution, nye bygninger"
    F         "Substitution, nye og gamle bygninger"
    fp          "Priskorrektion"
    fq          "Vækstkorrektion"
    fv          "Pris- og vækstkorrektion"
    beta      "Beta"
    rho         "Rho"
    r         "Realrente"
    c         "Løbende omkostninger"
;

Equations
E_1(a,t)
E_2(a,t)
E_3(t)           
E_4(a,t)        
E_5(t)             
E_6(t)           
E_7(t)                     
E_8(t)
*E_c_9
E_10
E_11
E_12
E_13
E_14
E_15
E_16
E_17
*E_18
;

*------------------- Modelligninger:

E_1(a,t-1)$axA(a)..      p_uc(a,t-1)                  =E= r(t)*p(a,t-1)/fp + (p(a,t-1)/fp - p(a+1,t)) + c(a,t);
*E_18(a,t-1)$aA(a)..      p_uc(a,t-1)                  =E= p_uc(a-1,t-1);


E_2(a,t-1)$axA(a)..      K(a,t-1)                     =E= mu(a)*(p_uc(a,t-1)/(p_new(t-1)/fp))**(-E)*K_new(t-1);
E_3(t-1)..               p_new(t-1) * K_new(t-1)      =E= sum(axA,p_uc(axA,t-1)*K(axA,t-1));
E_4(a,t-1)$aA(a)..       p_old(t-1)                   =E= (r(t)*p(a,t-1)/fp+(p(a,t-1)/fp-p(a,t))+c(a,t))*fp;

*Ligninger fra noten
E_5(t-1)$txT(t)..             K_new(t-1)    =E= gamma_new*(p_new(t-1)/p_U(t-1))**(-F)*u(t);
E_6(t-1)$txT(t)..             K_old(t-1)    =E= gamma_old*(p_old(t-1)/p_U(t-1))**(-F)*u(t);

E_16(t-1)$tT(t)..                    K_new(t) =E= K_new(t-1);
E_17(t-1)$tT(t)..                    K_old(t) =E= K_old(t-1);

E_7(t-1)..                      p_U(t) * u(t)        =E= (p_new(t-1)*K_new(t-1) + p_old(t-1)*K_old(t-1));
E_8(t-1)..                      W_hat(t) * fv            =E= (1+r(t))*W_hat(t-1) + Y_H(t) * fv - p_U(t-1)*u(t-1);
*E_c_9(t)..                               W_c_hat(t)                     =E= W_c(t) + sum(a, p_c(a,t)*K_c(a,t));

E_10(t-1)$txT(t)..              (u(t) / fq)**(-rho)     =E= (((p_U(t)/fp)/p_U(t+1))*(1+r(t+1))*Beta*u(t+1)**(-Rho));
E_11(t-1)$tT(t)..                                u(t)     =E= u(t-1);
E_12(a-1,t-1)$axA(a)..                          K(a,t)    =E= K(a-1,t-1)/fq;
E_13(a,t-1)$aA(a)..                           K_old(t-1)  =E= K_old(t-1)/fq+K(a-1,t-1)/fq-N_old(t);
E_14(a,t-1)$a0(a)..                           K(a,t-1)    =E= N_old(t-1) + G_bar(t-1);
E_15(a,t-1)$a0(a)..                           p(a,t-1)    =E= C_big(t-1) + p_old(t-1);

*------------------- Antagelser:
r.fx(t)        = 0.05;
c.fx(a,t)      = 1 + 0.01*ord(a);
F.fx           = 0.5;
E.fx           = 0.5;
fp.fx          = 1.00;
fq.fx          = 1.00;
fv.fx          = fp.l*fq.l;
rho.fx         = 0.5;
beta.fx        = sum(t1,1/(1+r.l(t1)));

*------------------- Initialisering:
* Priser
p.l(a,t)      = 1 - 0.01*ord(a); 
p_uc.l(a,t)   = sum(t0, r.l(t0)*p.l(a,t0)/fp.l + (p.l(a,t0)/fp.l - p.l(a+1,t0)) + c.l(a,t0));
p_new.l(t)    = 1;
p_old.l(t)    = 1;
p_U.l(t)      = 1;

* Bygninger
K.l(a,t)     = 1;
K_new.l(t)   = sum(axA,p_uc.l(axA,t)*K.l(axA,t)) / p_new.l(t);
K_old.l(t)   = 200;
N_old.l(t)   = 0;

*Andet
u.l(t)     = (p_new.l(t) * K_new.l(t) + p_old.l(t) * K_old.l(t))/(p_U.l(t));

display K.l, K_new.l, K_old.l, u.l, p.l, p_uc.l;

*------------------- Eksogene input
G_bar.fx(t) = 1;
C_big.fx(t) = 1;
Y_H.fx(t)   = 1;

*K_old.fx(t0) = K_old.l(t0);
*K_new.fx(t0) = K_new.l(t0);

*N_old.fx(t0) = 0;
*P_old.fx(t0) = 1;
K.fx(a,t0) = K.l(a,t0);
p.fx(a,t0) = p.l(a,t0);

*------------------- Kalibrering:
mu.fx(a)     = K.l("0","2025")/((p_uc.l(a,"2025")/(p_new.l("2025")/fp.l))**(-E.l)*K_new.l("2025"));
gamma_new.fx = K_new.l("2025") / ((p_new.l("2025")/p_U.l("2025"))**(-F.l)*(p_new.l("2025") * K_new.l("2025") + p_old.l("2025") * K_old.l("2025"))/(p_U.l("2025")));
gamma_old.fx = K_old.l("2025")  / ((p_old.l("2025")/p_U.l("2025"))**(-F.l)*(p_new.l("2025") * K_new.l("2025") + p_old.l("2025") * K_old.l("2025"))/(p_U.l("2025")));

*w_c_hat.fx(t)  = (p_c_new.l("2024") * K_c_new.l("2024") + p_c_old.l("2024") * K_c_old.l("2024") - Y_H.l("2024"))/((1+r_c("2024"))/fv - 1);
w_hat.fx(t0)  = (p_new.l(t0) * K_new.l(t0) + p_old.l(t0) * K_old.l(t0) - Y_H.l(t0))/((1+r.l(t0))/fv.l - 1);
*w_c.fx(t)      = (W_c_hat.l("2024") +  K_c.l("0","2025") * sum(axA.local, P_c.l(axA,"2024")) + P_c.l("1","2024") * K_c_old.l("2024"));
c_big.fx(t)  = P.l("0","2024") - p_old.l("2024");

*------------------- Modelløsning

*model lille_model /E_c_1, E_c_2, E_c_3, E_c_4, E_c_5, E_c_6, E_c_7/; 

model model_byg /ALL/;
solve model_byg using CNS;


*p_c.fx(a,t) = p_c.l(a,t);
*solve lille_model using CNS;

