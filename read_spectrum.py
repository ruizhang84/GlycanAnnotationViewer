import numpy as np
import pandas as pd
from sklearn.cluster import KMeans

class SPECTRUM:
    def __init__(self, peaks, charge, scanno, mz):
        self.peaks = peaks
        self.charge = charge
        self.scan = scanno
        self.mz = mz

def read_mgf(fname):
    """
        input mgf
        output: spectrum list
    """
    spectra = {}
    charge = 0
    premass = 0.0
    scan = 0
    peaks = []
    with open(fname, 'r') as f:
        for line in f:
            entry = line.strip()
            if "BEGIN IONS" in entry:
                peaks = []
            elif "MASS=mono" in entry or "TITLE" in entry or "RTINSECONDS" in entry:
                continue
            elif "PEPMASS" in entry:
                temp = entry.split("=")
                temp = temp[1].split(" ")
                premass = float(temp[0])
            elif "CHARGE" in entry:
                temp = entry.split("=")
                charge = int(temp[1][:-1])
            elif "SCAN" in entry:
                temp = entry.split("=")
                scan = int(temp[1])

            elif "END IONS" in entry:
                df_peaks = pd.DataFrame(peaks)
                spec = SPECTRUM(df_peaks, charge, scan, premass)
                spectra[scan] = spec
            else:
                temp = entry.split(" ")
                if len(temp) == 2:
                    mz, intensity = temp
                    peaks.append({"intensity": float(intensity), "mz": float(mz)})
    return spectra

def insert_peaks(peaks, tol = 0.1):
    new_peaks = []
    prev = None
    for index, row in peaks.iterrows():
        mz = row.mz
        intensity = row.intensity
        if prev is None:
            prev = mz
        elif abs(prev - mz) > 0.1:
            new_peaks.append({"intensity": 0, "mz": prev + tol/2.0})
            new_peaks.append({"intensity": 0, "mz": mz - tol/2.0})
            prev = mz
        new_peaks.append({"intensity": intensity, "mz": mz})
    return pd.DataFrame(new_peaks)        

def cluster_peaks(peaks, n_clusters):
    X = np.array(peaks.intensity).reshape(-1, 1)
    kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(X)
    
    low_index = 0
    low_intensity = float("inf")
    for i in range(n_clusters):
        val = peaks[kmeans.labels_ == i].intensity.mean()
        if val < low_intensity:
            low_intensity = val
            low_index = i
    return kmeans, low_index