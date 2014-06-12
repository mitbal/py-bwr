# Module DSP, Collection of method for digital signal processing

import math
import random

# Input Side Convolution
def conv(x, h):

    length = len(x) + len(h) - 1
    y = [0]*length

    for i in xrange(len(x)):
        for j in xrange(len(h)):
            y[i+j] += x[i] * h[j]

    return y

# Output Side Convolution
def conv2(x, h):

    length = len(x) + len(h) - 1
    y = [0]*length

    for i in xrange(len(y)):
        for j in xrange(len(h)):
            if i-j >= 0 and i-j < len(x):
                y[i] += h[j] * x[i-j]

    return y

# Discrete Fourier Transform
def dft(x) :
    N = len(x)
    X = [0]*N

    for k in xrange(N) :
        for n in xrange(N) :
            X[k] += x[n] * math.e**(-2j * math.pi * k * n / N)

    return X

# Fast Fourier Transform
def fft(x):
    N = len(x)
    if N == 1:
        return x

    x_odd = x[::2]
    x_even = x[1::2]

    X_odd = fft(x_odd)
    X_even = fft(x_even)

    # Butterfly
    X = [0]*N
    half = N/2
    for i in xrange(half):
        X[i] = X_even[i] + X_odd[i]*math.e**(-2j*math.pi*i/N)
        X[i+half] = X_even[i] - X_odd[i]*math.e**(-2j*math.pi*i/N)

    return X

# Fast Fourier Transform with real and imaginary part of the number separated
def fft2(x_re, x_im):

    N = len(x_re)
    if N == 1:
        return (x_re, x_im)

    x_even_re = x_re[::2]
    x_odd_re = x_re[1::2]

    x_even_im = x_im[::2]
    x_odd_im = x_im[1::2]

    (X_even_re, X_even_im) = fft2(x_even_re, x_even_im)
    (X_odd_re, X_odd_im) = fft2(x_odd_re, x_odd_im)

    # Butterfly
    X_re = [0]*N
    X_im = [0]*N
    half = N/2
    for k in xrange(half):
        tr = X_odd_re[k]*math.cos(-2*math.pi*k/N) - X_odd_im[k]*math.sin(-2*math.pi*k/N)
        ti = X_odd_im[k]*math.cos(-2*math.pi*k/N) + X_odd_re[k]*math.sin(-2*math.pi*k/N)

        X_re[k] = X_even_re[k] + tr
        X_im[k] = X_even_im[k] + ti

        X_re[k+half] = X_even_re[k] - tr
        X_im[k+half] = X_even_re[k] - ti

    return (X_re, X_im)

# Inverse Fast Fourier Transform
def ifft(X):

    N = len(X)
    for i in xrange(N):
        X[i] = -1j*X[i]

    x = fft(X)

    for i in xrange(N):
        x[i] = x[i]/N*-1j

    return x

# Polar coordinates of frequency domain
def polar(X):

    N = len(X)
    mag = [0]*N
    phase = [0]*N

    for i in xrange(N):
        mag[i] = abs(X[i])
        phase[i] = math.arctan(X[i].imag/X[i].real)

    return (mag, phase)

# Convolution 2 Dimension
def conv2d(x, h):
    height = len(x)
    width = len(x[0])

    y = [None]*height
    for i in xrange(len(y)):
        y[i] = [0]*width

    return y;


# Moving Average Filter Implemented with convolution
def mav(x, m):

    h = [1.0/m]*m
    y = conv(x, h)
    return y

# Moving Average Filter implemented with recursion
def mav2(x, m):
    y = [1.0/m]*m

    return y

# Method for generating range in floating point
def frange(start, end, num):
    x = [0]*num;
    diff = end-start
    step = diff/(num-1)

    for i in xrange(num):
        x[i] = start+(step*i)

    return x

# Generate sin function
def gen_sine(num, k):
    sinx = [0]*num
    for i in xrange(num):
        sinx[i] = math.sin(k*2*math.pi*(i/float(num)))

    return sinx

# Generate cos function
def gen_cosine(num, k):
    cosx = [0]*num
    for i in xrange(num):
        cosx[i] = math.cos(k*2*math.pi*(i/float(num)))

    return cosx

# Combine two signals using addition
def add(x, y):
    z = [0]*len(x)
    for i in xrange(len(x)):
        z[i] = x[i]+y[i]

    return z

# Combine two signals using multiplication
def mul(x, y):
    N = len(x)
    z = [0]*N
    for i in xrange(N):
        z[i] = x[i]*y[i]

    return z

# Generate random signal with maximum and minimum value
def gen_random(num, max, min):
    x = [0]*num
    scope = max-min

    for i in xrange(num):
        x[i] = random.random()*scope+min

    return x

# Generate low pass filter kernel
def gen_low_pass_kern(fc, m):
    h = [0]*m

    # Generate kernel filter
    for i in xrange(m):
        point = i-m/2
        if point == 0:
            h[i] = 2*math.pi*fc
        else:
            h[i] = math.sin(2*math.pi*fc*point)/point
        h[i] = h[i]*(0.42 - 0.5*math.cos(2*math.pi*i/m) + 0.08*math.cos(4*math.pi*i/m))

    # Normalization
    sum = 0
    for i in xrange(m):
        sum += h[i]
    for i in xrange(m):
        h[i] /= sum

    return h

# Generate high pass kernel
def gen_high_pass_kern(fc, m):

    # Generate low pass first
    hl = gen_low_pass_kern(fc, m)

    # Spectral Inversion
    for i in xrange(m):
        hl[i] = -hl[i]
    hl[m/2] = 1-hl[m/2]

    return hl

# Implementation of low pass filter windowed-sinc with Blackman Window
def low_pass(x, fc, m):

    # Generate low pass filter kernel
    h = gen_low_pass_kern(fc, m)

    y = conv(x, h)
    return y

# Implementation of High pass filter
def high_pass(x, fc, m):

    # Generate the high pass kernel
    h = gen_high_pass_kern(fc, m)

    # The actual convolution
    y = conv(x, h)
    return y

# Implementation of Band pass filter
def band_pass(x, fcmin, fcmax, m):

    # Creating the low pass filter
    hl = gen_low_pass_kern(fcmin, m)

    # Creating the high pass filter
    hh = gen_high_pass_kern(fcmax, m)

    # Convolve the two filter together
    h = conv(hl, hh)

    # Convolve with the input signal
    y = conv(x, h)
    return y

# Implementation of Band reject filter
def band_reject(x, fcmin, fcmax, m):

    # Creating the low pass filter
    hl = gen_low_pass_kern(fcmin, m)

    # Creating the high pass filter
    hh = gen_high_pass_kern(fcmax, m)

    # Combine the two kernels
    h = combine(hl, hh)

    # Convolve with the actual signal
    y = conv(x, h)
    return y

# Filter with recursive implementation
def rec_filter(x, a, b):
    N = len(x)
    y = [0]*N

    for i in xrange(N):
        y[i] = a[0]*x[i] + a[1]*x[i-1] + a[2]*x[i-2] + b[1]*y[i-1] +b[2]*y[i-2]

    return y

# FFT Convolution
def fft_conv(x, h):

    # FFT input signal
    X = fft(x)

    # FFT kernel filter
    H = fft(h)

    # Multiplication of signal
    N = len(X)
    Y = [0]*N
    for i in xrange(N):
        Y[i] = X[i]*H[i]

    # Inverse FFT
    y = ifft(Y)
    return y

# Convolution through overlap and add
def overlap_add(x, h):
    y = len(X)+len(h)-1

    return y
