# -*- coding: utf-8 -*-

import numpy as np

def FourierPowerSpectrum(acc: np.ndarray, dt: float):
    """
    acc  : ground acceleration time history, unit = g
    dt   : time step, unit = s

    return:
    freq  : frequency, unit = Hz
    amp   : Fourier amplitude spectrum, unit = cm/s
    power : power spectrum (square spectrum), unit = cm^2/s^2
    """

    n = len(acc)

    # ---------------------------------------------------
    # g → cm/s²
    # ---------------------------------------------------
    acc_cms2 = acc * 981.0

    # 去均值
    acc0 = acc_cms2 - np.mean(acc_cms2)

    # FFT
    fft_vals = np.fft.rfft(acc0)
    freq = np.fft.rfftfreq(n, d=dt)   # Hz

    # ---------------------------------------------------
    # 傅氏幅值谱（cm/s）
    # ---------------------------------------------------
    amp = np.abs(fft_vals) / n
    if n > 1:
        amp[1:-1] *= 2.0

    # ---------------------------------------------------
    # 功率谱（平方谱）（cm²/s²）
    # ---------------------------------------------------
    power = amp ** 2

    return freq, amp, power


if __name__ == "__main__":
    pass