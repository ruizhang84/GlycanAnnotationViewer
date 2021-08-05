from tkinter import *
from tkinter import ttk 
from drawer import Drawer
from PIL import ImageTk, Image

class MakeWindow:
    def __init__(self, root,
        f1="10MixGlycanStandards_C18_50cm_091520.mgf", 
        f2="10MixGlycanStandards_C18_50cm_091520_annotated.csv"):

        self.window = Toplevel(root) #instead of super
        self.window.title("Spectrum Annotations")
        self.window.geometry("1500x600")

        self.drawer = Drawer(f1, f2)
        print (f1, f2)
        self.scan_select = IntVar()
        self.display_info = StringVar()
        self.display_info.set("scan: ")

        label = ttk.Label(self.window, textvariable=self.display_info).grid(row=0, column=0, padx=5, pady=10)
        scale = ttk.Scale(self.window, command=lambda x: self.call_back(x, self.drawer.scans), orient= HORIZONTAL, length= 700, 
            variable=self.scan_select, from_ = 0, to = len(self.drawer.scans)-1).grid(row=0, column=1, padx=5, pady=10)

        button_annot = ttk.Button(self.window, command=lambda :self.plot(self.drawer), text="View").grid(row=0, column=2, padx=5, pady=10)

    def call_back(self, x, scans):
        index = int(float(x))
        info = f"scan: {scans[index]}"
        self.display_info.set(info)

    def plot(self, drawer):
        index = int(self.scan_select.get())
        scan = self.drawer.scans[index]

        # plot glycan spectrum image
        fig = self.drawer.draw(scan, figsize=(14, 5))

        # plot glycan structure
        img = Image.open("glycan.png")
        (width, height) = (img.width // 3, img.height // 3)
        img = img.resize((width, height), Image.ANTIALIAS)
        img = ImageTk.PhotoImage(img)
        panel = Label(self.window, image=img)
        panel.image = img
        panel.grid(row=1, column=0, pady=10)

        # plot spectrum
        img = Image.open("annot.png")
        img = ImageTk.PhotoImage(img)
        panel = Label(self.window, image=img)
        panel.image = img
        panel.grid(row=1, column=1, rowspan=4, columnspan=4, pady=10)

