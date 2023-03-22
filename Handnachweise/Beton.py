#Verwendete Packages
from math import *
from handcalcs import handcalc

import forallpeople
forallpeople.environment('structural', top_level=True)


#Nachweise für Platten

## Bestimmung des Rissmoments

@handcalc('short',jupyter_display=True, precision=2)
def platte_rissmoment(h, f_ctm):
    t = h/3
    k_t = 1/(1+0.5*t.value) #h in Meter
    f_ctd = k_t*f_ctm
    m_cr = h**2/6*f_ctd
    return locals()

## Bestimmung des Biegewiderstands einer Platte

@handcalc('long', jupyter_display=True, precision=2)
def platte_biegewiderstand(D, s, h, c_nom, f_sd, f_cd):
    a_s = pi*(D/2)**2/s
    d_s = h-c_nom-D-D/2
    m_Rd = a_s*f_sd * (d_s-(a_s*f_sd)/(2*f_cd))
    x =  a_s * f_sd /(0.85*f_cd)

    Kontrolle = x/d_s <=0.35
    return locals()

## Betonquerkraftwiderstand

@handcalc('long', jupyter_display=True, precision=2)
def platte_betonquerkraft(k_c, f_cd, d_v, theta_c3):

    v_Rdc = k_c * f_cd * d_v * sin(theta_c3)*cos(theta_c3)
    return locals()

## Plattenquerkraft mit Bewehrung

@handcalc(jupyter_display=True, precision=3)
def platte_querkraft_normativ(v_xd, v_yd, m_xd, m_yd, m_xyd, m_xRd, m_yRd, f_sd, E_s, d_v, D_max, tau_cd):
    v_0d = sqrt(v_xd**2 + v_yd**2)
    phi_0d = atan(v_yd/v_xd)
    phi_0_d_deg = (degrees(phi_0d))

    m_nRd = m_xRd * cos(phi_0d)**2 + m_yRd*sin(phi_0d)**2
    m_nd = cos(phi_0d)**2*m_xd + sin(phi_0d)**2*m_yd + sin(2*phi_0d)*m_xyd

    k_g = 48/(16+D_max)* m_nd/m_nd

    zeta = 1/(sin(phi_0d)**4+cos(phi_0d)**4)

    e_v_EL = f_sd/E_s*abs(m_nd/m_nRd) 
    e_v_ZL = 1.5*f_sd/E_s 

    k_d_EL = 1/(1+e_v_EL*float(d_v)*k_g*zeta)
    k_d_ZL = 1/(1+e_v_ZL*float(d_v)*k_g*zeta)

    v_Rd_EL = k_d_EL * tau_cd * d_v 
    v_Rd_ZL = k_d_ZL * tau_cd * d_v
    return locals()

@handcalc(jupyter_display = True, precision=3,override='long')
def platte_querkraftwiderstand_stuetze(h,tau_cd,f_cd,d_max,E_s,f_sd,c_nom,c_s,phi_x,phi_y,a_x,a_y,a,b,m_sdx,m_sdy,r_sx,r_sy,k_e):
    k_g = (48) / (16+d_max) # 37
        
    d_x = h-c_nom-phi_x/2
    d_y = h-c_nom-phi_x-phi_y/2

    d_v = (d_x + d_y)/2-c_s # Fig. 20
        
    a_sx = (phi_x/2)**2*pi*1000/a_x #Bewehrungsquerschnitt in x-Richtung
    a_sy = (phi_y/2)**2*pi*1000/a_y #Bewehrungsquerschnitt in y-Richtung
            
    m_Rdx = a_sx*f_sd*(d_x-a_sx*f_sd/(2*1000*f_cd))
    m_Rdy = a_sy*f_sd*(d_y-a_sy*f_sd/(2*1000*f_cd))
           
    psi_x = 1.5 * r_sx / d_x * f_sd / E_s * (m_sdx/m_Rdx)**(3/2) # 59
    psi_y = 1.5 * r_sy / d_y * f_sd / E_s * (m_sdy/m_Rdy)**(3/2) # 59
        
    k_rx = 1 / (0.45 + 0.18 * psi_x * d_x * k_g) # 58
    k_ry = 1 / (0.45 + 0.18 * psi_y * d_y * k_g) # 58
        
    k_r = min(k_rx,k_ry) # Kleinerer Widerstand der beiden Richtungen

    u_e = 2*(a+b)+d_v*pi # Umfang des Nachweisschnitts
        
    V_Rdc = k_r*tau_cd*d_v*u_e*k_e
    return locals()
