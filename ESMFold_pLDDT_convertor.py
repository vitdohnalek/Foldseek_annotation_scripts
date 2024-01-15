#####################################################
#    Changes the pLDDT values in bfactor column     #
#    Needed for the EMSFold files from esmatlas     #
#    Before you color them in ChimeraX              #
#####################################################

import tkinter as tk
from tkinterdnd2 import TkinterDnD, DND_FILES
import os

def process_pdb_file(file_path):
    new_pdb = ""
    
    with open(file_path, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                plddt = int(float(line[62:66]) * 100)
                new_pdb += line[:62] + str(plddt) + line[66:]
            else:
                new_pdb += line
    
    new_file_path = f"{file_path[:-4]}_conv.pdb"
    with open(new_file_path, "w") as f:
        f.write(new_pdb)
    
    return new_file_path

def handle_drop(event):
    file_path = event.data
    if file_path.lower().endswith('.pdb'):
        new_file_path = process_pdb_file(file_path)
        result_label.config(text=f"Updated PDB file created:\n {new_file_path}")

# Create the main window
root = TkinterDnD.Tk()

#Set size of the window
root.geometry("600x150")

# Set the title of the window
root.title("ESMFold pLDDT Converter")

# Set up the drag-and-drop functionality
root.drop_target_register(DND_FILES)
root.dnd_bind('<<Drop>>', handle_drop)

# Create a label for displaying the result
result_label = tk.Label(root, text="Drag and drop a .pdb file here")
result_label.pack(pady=50)

# Start the main loop
root.mainloop()
