from mg5.CalcMEZHH import lib as zhh
from mg5.CalcMEZZH import lib as zzh
from Physsim.CalculateME import lib as Physsim

fm = [
    20, 1, 2, 3,
    25, 2, 3, 4,
    100, 3, 4, 5,
    120, 4, 5, 6,
    150, 5, 6, 7,
    160, 6, 7, 8, 
]

fm_ph = [
    20, 1, 2, 3,
    25, 2, 3, 4,
    220, 7, 9, 11,
    310, 11, 13, 15
]

print(zhh.calc_single(b"/afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/analysis/cffi/mg5/mg5/Cards/param_card.dat", 2, [0,-1,1,1], fm))
print(Physsim.calc_zhh_single(fm_ph))

#zhh.calc(b"/afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/analysis/cffi/mg5/mg5/Cards/param_card.dat", 2, [0,-1,1,1])
#zzh.calc(b"/afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/analysis/cffi/mg5/mg5/Cards/param_card.dat", 2, [0,-1,1,1])