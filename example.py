import bwr
import matplotlib.pyplot as plt

# Read input csv file from physionet
f = open('samples1.csv', 'r')
lines = f.readlines()
f.close()

# Discard the first two lines because of header. Takes either column 1 or 2 from each lines (different signal lead)
signal = [0]*(len(lines)-2)
for i in xrange(len(signal)):
	signal[i] = float(lines[i+2].split(',')[1])

# Call the BWR method
(baseline, ecg_out) = bwr.bwr(signal)

plt.subplot(2,1,1)
plt.plot(signal, 'b-')
plt.plot(baseline, 'r-')

plt.subplot(2,1,2)
plt.plot(ecg_out, 'b-')
plt.show()
