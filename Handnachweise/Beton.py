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