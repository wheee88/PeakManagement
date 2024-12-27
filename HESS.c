#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h> // POSIX 표준 sleep 함수를 위해

#define INF 1000000000.0 // 양의 무한대
double OnSig = 0;
double dmode = 0;
double hmode = 0;
double cmode = 0;
double mode = 0;
double TInt = 0.001;
double SC_UL = 0.0;
double SC_LL = 0.0;
double BE_UL = 0.0;
double BE_LL = 0.0;
double BESS_ramp_lim = 0.0;
double P_RefB = 0.0;
double P_ACLoad = 0.0;
double P_mng = 0.0;
double Ptot;
double E_sub;
double E_manage;
double E_diff;
double P_err_SCES;
double P_Int_SCES;
double P_pi_SCES = 1.65;
double P_err_BESS;
double P_Int_BESS;
double P_pi_BESS = 1.65;

double Ki_P_SCES = 0.01 / 50e-6;
double Kp_P_SCES = 0.01;
double Ki_P_BESS = 0.1 / 50e-6;
double Kp_P_BESS = 0.001;
static double prev_Vsrc_BESS = 0.0;
static double ntimes = 0.0;
static double reset = 0.0;
static double PD_SCES;
static double PH_SCES;
static double PC_SCES;
static double PD_BESS;
static double PH_BESS;
static double PC_BESS;
static int InitFlag = 0;

double NOR(double A, double B) {
    return !((A != 0) || (B != 0)); // A와 B가 모두 0인 경우에만 1 반환
}

double Limiter(double x, double U_limit, double L_limit) // 열차제어 제한용
{
    double y = (x >= U_limit) ? U_limit : ((x <= L_limit) ? L_limit : x);
    return y;
}

// Rate Limiter 함수
double rate_limiter(double input, double prev_output, double up_rate, double down_rate, double delta_t) {
    double max_increase = up_rate * delta_t;   // 허용되는 최대 상승량
    double max_decrease = down_rate * delta_t; // 허용되는 최대 하강량
    double output;

    if (input > prev_output + max_increase) {
        output = prev_output + max_increase; // 상승 제한
    } else if (input < prev_output - max_decrease) {
        output = prev_output - max_decrease; // 하강 제한
    } else {
        output = input; // 제한 없이 그대로 출력
    }

    return output;
}

void hessctrl_(
    double *V_ESS_BUS, double *V_Load, double *V_char, double *V_disc, 
    double *P_SUB, double *P_Rail, double *P_ACLoad, double *P_BESS, 
    double *P_SCES, double *E_BESS, double *E_SCES, double *R_Sub, double *R_int, double *Tlength, 
    double *BESS_Capa, double *BESS_ramp, double *SCES_Capa, double *Ppeak, 
    double *managed, double *Vsrc_SCES, double *Vsrc_BESS, double *Timer, 
    double *Tdelt, double *Debug
) {
    if (*Timer >= (*Tlength * ntimes + 1.0)) 
    {
        if (reset != 1.0) { // reset이 아직 설정되지 않은 경우에만 실행
            reset = 1.0;
            ntimes += 1.0; // ntimes 증가
        }
    } else 
    {
        reset = 0.0; // 조건을 충족하지 않으면 reset을 초기화
    }

    // Operation delay switching signal to remove initial transient phenomenon
    if (*Timer >= 1.0) 
    {
        OnSig = 1.0;
    }

    // Operation mode selector
    dmode = (*V_disc >= *V_ESS_BUS) ? 1.0 : 0.0;
    dmode = dmode * OnSig;

    cmode = (*V_ESS_BUS >= *V_char) ? 1.0 : 0.0;
    cmode = cmode * OnSig;

    hmode = NOR(dmode, cmode);

    // Limiter boundary
    SC_UL = ((*SCES_Capa / *V_char) * *R_int) + *V_char;
    SC_LL = -((*SCES_Capa / *V_disc) * *R_int) + *V_disc;
    BE_UL = ((*BESS_Capa / *V_char) * *R_int) + *V_char;
    BE_LL = -((*BESS_Capa / *V_disc) * *R_int) + *V_disc;

    BESS_ramp_lim = (OnSig == 1) ? ((*BESS_ramp / *V_ESS_BUS) * *R_int) : ((*BESS_ramp / 1.65) * *R_int);

    // Substation capacity limit
    P_RefB = ((*V_char - *V_disc) / *R_Sub) * *V_disc;

    // Power reference for peak power management
    double temp_PSUB = *P_SUB;
    *P_SUB = Limiter(temp_PSUB, INF, 0.0);

    Ptot = (*P_ACLoad + *P_SUB);
    Ptot /= 3600.0;
    if (reset != 1.0) {
        E_sub += (TInt * Ptot);
    } else {
        E_sub = 0;
    }

    *Ppeak /= 3600.0;
    if (reset != 1.0) {
        E_manage += (TInt * *Ppeak);
    } else {
        E_manage = 0;
    }

    E_diff = Limiter((E_sub - E_manage), INF, 0.0) * 3.6;
    E_diff /= (*Tlength - *Timer + 1.0);

    P_mng = (E_diff * *managed * OnSig);

    if (!InitFlag)  // 초기 PI 제어 출력 설정
    {
        *Vsrc_SCES = 1.65;
        *Vsrc_BESS = 1.65;
        P_Int_SCES = 1.65;
        P_Int_BESS = 1.65;
        InitFlag = 1;
    }
    else 
    {
        PD_SCES = ((*P_SUB + P_mng - P_RefB) * dmode); // Discharge SCES
        PH_SCES = ((-*P_SCES + P_mng) * hmode); // Hold SCES
        PC_SCES = (((-*V_ESS_BUS + *V_char) / *R_Sub) * *V_ESS_BUS) * cmode; // Charge SCES

        P_err_SCES = PD_SCES + PH_SCES + PC_SCES;
        P_Int_SCES += (Ki_P_SCES * P_err_SCES * *Tdelt);
        P_pi_SCES = (Kp_P_SCES * P_err_SCES) + P_Int_SCES;

        PD_BESS = (-*P_BESS + *P_Rail - P_RefB) * dmode; // Discharge BESS
        PH_BESS = ((*managed == 1) ? (*P_SCES) : (-*P_BESS) ) * hmode; // Hold BESS
        PC_BESS = (*P_SCES + ((- *V_ESS_BUS + *V_char) / *R_Sub) * *V_ESS_BUS) * cmode; // Charge BESS

        P_err_BESS = PD_BESS + PH_BESS + PC_BESS;
        P_Int_BESS += (Ki_P_BESS * P_err_BESS * *Tdelt);
        P_pi_BESS = (Kp_P_BESS * P_err_BESS) + P_Int_BESS;
        *Vsrc_SCES = Limiter(P_pi_SCES, SC_UL, SC_LL);
        *Vsrc_BESS = Limiter(P_pi_BESS, BE_UL, BE_LL);
        *Vsrc_BESS = rate_limiter(*Vsrc_BESS, prev_Vsrc_BESS, BESS_ramp_lim, BESS_ramp_lim, *Tdelt);

    }

    prev_Vsrc_BESS = *Vsrc_BESS;

    if (reset != 1.0) {
        *E_SCES += (TInt * *P_SCES);
        *E_BESS += (TInt * *P_BESS);
    } else {
        *E_SCES = 0;
        *E_BESS = 0;
    }

    Debug[0] = P_pi_SCES;
    Debug[1] = P_pi_BESS;
    Debug[2] = PC_BESS;
    Debug[3] = *V_ESS_BUS;
    Debug[4] = SC_UL;
    Debug[5] = BESS_ramp_lim;
}
