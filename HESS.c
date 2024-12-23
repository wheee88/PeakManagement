#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h> // POSIX 표준 sleep 함수를 위해

int OnSig = 0;
int dmode = 0;
int hmode = 0;
int cmode = 0;
int mode = 0;
double SC_UL = 0.0;
double SC_LL = 0.0;
double BE_UL = 0.0;
double BE_LL = 0.0;
double BESS_ramp_lim = 0.0;
double P_RefB = 0.0;
double P_ACLoad = 0.0;
double P_mng = 0.0;
double P_Sub = 0.0;
double Vsrc_SCES = 0.0;
double Vsrc_BESS = 0.0;
double reset = 0.0;

// NOR 함수 정의
int NOR(int A, int B) {
    return !(A || B); // OR 결과의 NOT
}

void HESSCtrl_(double *V_ESS_BUS, double *V_Load,double *V_char,double *V_disc,double *P_SUB, double *P_Rail,double *P_Load,double *P_BESS,double *P_SCES,double *R_Sub,double *R_int,double *Tlength,double *BESS_Capa,double *BESS_ramp,double *SCES_Capa,double *Ppeak,double *CB_SC,double *CB_Bat,double *managed, double *Time, double *Debug)
{
    // Operation delay switching signal to remove initial transient phenomenon
    if(*Time > 1.0)
    {
        OnSig = 1;
    }

    // Operation mode selector
    if (*V_disc >= *V_ESS_BUS)
    {
        dmode = 1;
        dmode = dmode*OnSig;
    }
    if (*V_ESS_BUS >= *V_char)
    {
        cmode = 1;
        cmode = cmode*OnSig;
    }

    hmode = NOR(dmode,cmode);

    SC_UL = (*SCES_Capa / *V_char) * *R_int + *V_char; // Limiter boundary
    SC_LL = (*SCES_Capa / *V_disc) * *R_int + *V_disc; // Limiter boundary
    BE_UL = (*BESS_Capa / *V_char) * *R_int + *V_char; // Limiter boundary
    BE_LL = (*BESS_Capa / *V_disc) * *R_int + *V_disc; // Limiter boundary

    if (OnSig == 1)
    {
        BESS_ramp_lim = (*BESS_ramp / *V_ESS_BUS) * *R_int;
    }
    else
    {
        BESS_ramp_lim = (*BESS_ramp / 1.650) * *R_int;
    }

    P_RefB = ((*V_char - *V_disc) / *R_Sub) * *V_disc; // Substation capacity limit
}
