from math import sqrt
from dsp import conv

# Daubechies 4 Constant
c0 = (1+sqrt(3))/(4*sqrt(2))
c1 = (3+sqrt(3))/(4*sqrt(2))
c2 = (3-sqrt(3))/(4*sqrt(2))
c3 = (1-sqrt(3))/(4*sqrt(2))

def db4_dec(x, level):
    """
    Fungsi db4_rec melakukan dekomposisi wavelet dengan basis daubechies order 4 sebanyak level yang ditentukan
    """

    result = [[]]*(level+1)
    lpk = [c0, c1, c2, c3]
    hpk = [c3, -c2, c1, -c0]

    x_temp = x[:]
    for i in xrange(level):
        lp = conv(x_temp, lpk)
        hp = conv(x_temp, hpk)

        lp_ds=[0]*(len(lp)/2)
        hp_ds=[0]*(len(hp)/2)

        # Downsampling hingga menjadi setengah
        for j in xrange(len(lp_ds)):
            lp_ds[j] = lp[2*j+1]
            hp_ds[j] = hp[2*j+1]

        result[level-i] = hp_ds
        x_temp = lp_ds[:]

    result[0] = lp_ds
    return result

def db4_rec(signals, level):

    cp_sig = signals[:]
    # Reconstruction coefficient
    lpk = [c3, c2, c1, c0]
    hpk = [-c0, c1, -c2, c3]

    for i in xrange(level):
        lp = cp_sig[0]
        hp = cp_sig[1]

        # Verify new length
        length = 0
        if len(lp) > len(hp):
            length = 2*len(hp)
        else:
            length = 2*len(lp)

        # Upsampling
        lpu = [0]*(length+1)
        hpu = [0]*(length+1)
        index = 0
        for j in xrange(length+1):
            if j%2 != 0:
                lpu[j] = lp[index]
                hpu[j] = hp[index]
                index += 1

        # Convolution
        lpc = conv(lpu, lpk)
        hpc = conv(hpu, hpk)

        # Truncate
        lpt = lpc[3:-3]
        hpt = hpc[3:-3]

        # Add both signals
        org = [0]*len(lpt)
        for j in xrange(len(org)):
            org[j] = lpt[j] + hpt[j]

        if len(cp_sig) > 2:
            cp_sig = [org]+cp_sig[2:]
        else:
            cp_sig = [org]

    return cp_sig[0]

def norm(x):
    total = 0
    for i in x:
        total += i*i
    return sqrt(total)

def bwr(raw):

    en0 = 0
    en1 = 0
    en2 = 0
    n = 0

    curlp = raw[:]
    num_dec = 0
    last_lp = []
    while True:
        print 'Iterasi ke' + str(num_dec+1)
        print len(curlp)

        # Dekomposisi 1 level
        [lp, hp] = db4_dec(curlp,1)

        # Hitung energi detail
        en0 = en1
        en1 = en2
        en2 = norm(hp)
        print en2

        if en0 > en1 and en1 < en2:
            last_lp = curlp
            break

        curlp = lp[:]
        num_dec = num_dec+1

    # Hitung base-nya dengan rekonstruksi
    base = last_lp[:]
    for i in xrange(num_dec):
        base = db4_rec([base,[0]*len(base)], 1)

    # Kurangi sinyal asli dengan base-nya
    ecg_out = [0]*len(raw)
    for i in xrange(len(raw)):
        ecg_out[i] =  raw[i] - base[i]

    return (base, ecg_out)


