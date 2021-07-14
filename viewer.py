from tkinter import *
from tkinter import ttk 
from drawer import Drawer
from PIL import ImageTk, Image

window = Tk()
window.title("Spectrum Annotation")

def call_back(x, scans):
    index = int(float(x))
    info = f"scan: {scans[index]}"
    display_info.set(info)

def plot(drawer):
    index = int(scan_select.get())
    scan = drawer.scans[index]
    # plot glycan
    fig = drawer.display(scan)

    img = Image.open("glycan_temp.png")
    img = ImageTk.PhotoImage(img)
    panel = Label(window, image=img)
    panel.image = img
    panel.grid(row=1)

    # plot spectrum
    fig = drawer.draw(scan)
    img = Image.open("temp.png")
    img = ImageTk.PhotoImage(img)
    panel = Label(window, image=img)
    panel.image = img
    panel.grid(row=2)


f1 = "10MixGlycanStandards_C18_50cm_091520.mgf"
f2 = "10MixGlycanStandards_C18_50cm_091520_annotated.csv"

drawer = Drawer(f1, f2)
scan_select = IntVar()
display_info = StringVar()
display_info.set("scan: ")

label = ttk.Label(window, textvariable=display_info).grid(row=0, column=0, padx=5, pady=10)
scale = ttk.Scale(window, command=lambda x: call_back(x, drawer.scans), orient= HORIZONTAL, length= 400, variable=scan_select,
    from_ = 0, to = len(drawer.scans)-1).grid(row=0, column=1, rowspan=2, padx=5, pady=10)

button_annot = ttk.Button(window, command=lambda :plot(drawer), text="View").grid(row=0, column=3, padx=5, pady=10)

window.mainloop()