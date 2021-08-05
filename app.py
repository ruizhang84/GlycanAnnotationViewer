from tkinter import *
from tkinter import ttk 
from tkinter import filedialog
from viewer import MakeWindow

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

def call_back():
    MakeWindow(root, mgf_path_var.get(), csv_path_var.get())

label_ms = ttk.Label(root, text="MS/MS Spectra File").grid(row=0, column=0, padx=5)
button_ms = ttk.Button(root, command=open_mgf_file, text="Read MS/MS Files (*.mgf)").grid(row=0, column=1, padx=5, pady=10)

label_annot = ttk.Label(root, text="Annotation CSV File").grid(row=1, column=0, padx=5)
button_annot = ttk.Button(root, command=open_csv_file, text="Read Annotation Files (*.csv)").grid(row=1, column=1, padx=5, pady=10)

button_submit = ttk.Button(root, command=call_back, text="View").grid(row=2, column=2, padx=5, pady=10)

root.mainloop()
