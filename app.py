from tkinter import *
from tkinter import ttk 
from tkinter import filedialog

root = Tk()
root.title("GlycanViewer")
mgf_path_var = StringVar()
csv_path_var = StringVar()

def open_mgf_file():
    tf = filedialog.askopenfilename(
            title="Open MGF file", 
            filetypes=(("MGF File", "*.mgf"),)
            )
    mgf_path_var.set(tf)

def open_csv_file():
    tf = filedialog.askopenfilename(
            title="Open CSV file", 
            filetypes=(("CSV File", "*.csv"),)
            )
    csv_path_var.set(tf)

def call_back2(x):
    print (x)

def call_back():
    window = Toplevel(root)
    window.title("Spectrum Annotations")
    
    scale = ttk.Scale(window, command=call_back2, orient= HORIZONTAL, length= 400, from_ = 1, to = 11).pack()
    print (mgf_path_var.get())
    print (csv_path_var.get())

class MakeWindow(tk.Toplevel):
    def __init__(self, message):
        tk.Toplevel.__init__(self) #instead of super
        self.message = message
        self.display = tk.Label(self, text=message)
        self.display.pack()

label_ms = ttk.Label(root, text="MS/MS Spectra File").grid(row=0, column=0, padx=5)
button_ms = ttk.Button(root, command=open_mgf_file, text="Read MS/MS Files (*.mgf)").grid(row=0, column=1, padx=5, pady=10)

label_annot = ttk.Label(root, text="Annotation CSV File").grid(row=1, column=0, padx=5)
button_annot = ttk.Button(root, command=open_csv_file, text="Read Annotation Files (*.csv)").grid(row=1, column=1, padx=5, pady=10)

button_submit = ttk.Button(root, command=call_back, text="View").grid(row=2, column=2, padx=5, pady=10)

root.mainloop()