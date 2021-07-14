from tkinter import *
from tkinter import ttk 

window = Tk()
window.title("Spectrum Annotation")

def call_back(x):
    print (x)


f1 = "10MixGlycanStandards_C18_50cm_091520.mgf"
f2 = "10MixGlycanStandards_C18_50cm_091520_annotated.csv"


scale = ttk.Scale(window, command=call_back, orient= HORIZONTAL, 
    length= 400, from_ = 1, to = 11).grid(row=0, column=1, rowspan=2, padx=5, pady=10)

window.mainloop()