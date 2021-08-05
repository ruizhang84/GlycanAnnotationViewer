import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

from matplotlib.patches import Circle
from matplotlib.offsetbox import (TextArea, DrawingArea, OffsetImage,
                                  AnnotationBbox)
from plot_glycan import glycan_table_plot
import read_spectrum
from read_spectrum import read_mgf, insert_peaks, cluster_peaks


def get_figure_size(glycan_id, size=4, line_length = 1000):
    fig, ax = plt.subplots(figsize=(size, size))
    glycan_table_plot(0, 0, line_length, ax, glycan_id)
    ax.set_axis_off()
    x1, x2 = ax.get_xlim()
    y1, y2 = ax.get_ylim()
    ratio = abs(x1 - x2) / abs(y1 - y2)
    fig.clear()
    ax.clear()
    plt.close()
    return  abs(x1 - x2) / line_length / 4, abs(y1 - y2) / line_length / 4

def generate_glycan(glycan_id, img_name="glycan.png", size=4, line_length = 1000):
    fig, ax = plt.subplots(figsize=get_figure_size(glycan_id, size, line_length))
    glycan_table_plot(0, 0, 1000, ax, glycan_id)
    ax.set_axis_off()
    fig.savefig(img_name, transparent=True)
    fig.clear()
    ax.clear()
    plt.close()

def get_glycan_id(glycan_string):
    return [int(i) for i in glycan_string.split()]

def get_fragment_id(fragment_string):
    glycan_id_string = fragment_string.split("|")[0].split(":")[-1]
    return get_glycan_id(glycan_id_string)

class Drawer:
    def __init__(self, 
        spectrum_path="10MixGlycanStandards_C18_50cm_091520.mgf", 
        annotated_path="10MixGlycanStandards_C18_50cm_091520_annotated.csv"):
        self.spectra = read_mgf(spectrum_path)
        self.df_mark = pd.read_csv(annotated_path)
        self.scans = sorted(self.df_mark.scan)

    def draw(self, scan=15666, figsize=(16, 6)):
        df_select = self.df_mark[self.df_mark.scan == scan]
        glycan_id = get_glycan_id(df_select.glycan.iloc[0])

        fig, ax = plt.subplots(figsize=figsize)

        kmeans, low_index = cluster_peaks(self.spectra[scan].peaks, 3)
        peaks = insert_peaks(self.spectra[scan].peaks)

        # plot spectrum
        ax.plot(np.round(peaks.mz, 2), peaks.intensity, "k")
        ax.set_xlabel("M/Z")
        ax.set_ylabel("Intensity")
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        # plot glycan
        generate_glycan(glycan_id)
        # img_x = np.min(peaks.mz) * 1.2
        # img_y = np.max(peaks.intensity)
        # img_insert = plt.imread('glycan.png')
        # imagebox = OffsetImage(img_insert, zoom=0.3)
        # imagebox.image.axes = ax

        # ab = AnnotationBbox(imagebox, (img_x, img_y))
        # ab.patch.set_edgecolor('none')
        # ab.patch.set_facecolor('none')
        # ax.add_artist(ab)

        # annotations
        for index, row in self.df_mark[self.df_mark.scan == scan].iterrows():
            mz = np.round(row.mz, 2)
            intensity = row.intensity
            xy = (mz, intensity)
            ax.plot([mz, mz], [0, intensity], "r")
            
            if kmeans.predict([[intensity]])[0] == low_index:
                continue
            
            generate_glycan(get_fragment_id(row.fragments), 'annot.png')
            img_insert = plt.imread('annot.png')
            imagebox = OffsetImage(img_insert, zoom=0.25)
            imagebox.image.axes = ax

            ab = AnnotationBbox(imagebox, (mz, intensity),
                            xybox=(22, np.log2(intensity)),
                            xycoords='data',
                            boxcoords="offset points",
                            arrowprops=dict(arrowstyle="->")
                            )
            ab.patch.set_edgecolor('none')
            ab.patch.set_facecolor('none')
            ax.add_artist(ab)
            
            offsetbox = TextArea(str(mz))

            ab = AnnotationBbox(offsetbox, (mz, intensity),
                            xybox=(-22, 2 * np.log2(intensity)),
                            xycoords='data',
                            boxcoords="offset points",
                            arrowprops=dict(arrowstyle="->")
                            )
            ab.patch.set_edgecolor('none')
            ab.patch.set_facecolor('none')
            ax.add_artist(ab)
        # save
        fig.savefig("annot.png", transparent=True)
        return fig

if __name__ == "__main__":
    draw = Drawer()
    draw.draw()