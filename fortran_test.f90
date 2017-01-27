Program fortran_test
a = 0.5_16
w_lam = -0.5_16
wa_ppf = 0.5_16
f_de = 0._16
f_de=a**(1.-3.*w_lam-3.*wa_ppf)*&
        exp(-wa_ppf*(1.-a)*(60.*a**6-430.*a**5+1334.&
        *a**4-2341.*a**3+2559.*a**2-1851.*a+1089.)/140.)
print *, f_de
END Program