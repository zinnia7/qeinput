#!/usr/bin/python3

import os
from ase.io import read
from ase.visualize import view
from ase.constraints import FixAtoms
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np

sections = ["&CONTROL","&SYSTEM","&ELECTRONS","&IONS","&CELL","CELL_PARAMETERS","ATOMIC_SPECIES","ATOMIC_POSITIONS","K_POINTS"]
keyword_defaults = {
    "&CONTROL": {
        "calculation": ["'scf'","'relax'","'vc-relax'","'nscf'","'bands'"],
        "prefix": "'pwscf'",
        "etot_conv_thr": 1e-4,
        "forc_conv_thr": 1e-3,
        "tefield": [".FALSE.",".TRUE."],
        "dipfield": [".FALSE.",".TRUE."]
    },
    "&SYSTEM": {
        "ibrav": 0,
        "nat": 0,
        "ntyp": 0,
        "ecutwfc": 30,
        "ecutrho": 180,
        "occupations": ["'smearing'","'tetrahedra'","'fixed'"],
        "smearing": ["'gauss'","'mp'","'mv'","'fd'"],
        "degauss": 0.01,
        "nspin": [1,2],
        "tot_magnetization": 0,
        "starting_magnetization(1)": [-1,1],
        "edir": [1,2,3],
        "emaxpos": 0.5,
        "eopreg": 0.1
    },
    "&ELECTRONS": {
        "electron_maxstep": 100,
        "conv_thr": 1e-6,
        "mixing_beta": 0.7,
        "mixing_mode": ["'plain'","'TF'","'local-TF'"],
        "diagonalization": ["'david'","'cg'","'ppcg'","'paro'","'rmm-davidson'"],
        "startingpot": ["'atomic'","'file'"],
        "startingwfc": ["'atomic'","'atomic+random'","'random'","'file'"]
    },
    "&IONS": {
        "ion_dynamics": ["'bfgs'","'damp'"]
    },
    "&CELL": {
        "cell_dynamics": ["'bfgs'","'damp-w'","'damp-pr'"],
        "press_conv_thr": 0.5,
        "cell_dofree": ["'all'","'shape'","'2Dxy'"]
    },
    "K_POINTS": {
        "K_POINTS": ["automatic 1 1 1 0 0 0", "gamma"],
    }
}

def load_input_file():
    file_path = filedialog.askopenfilename(filetypes=[("Input files", "*.pwi *.in")])
    if not file_path:
        return

    with open(file_path, "r") as input_file:
        lines = input_file.readlines()

    atoms = read(file_path, format="espresso-in")

    current_section = None
    atomic_positions = []
    ase_constraints = []
    for line in lines:
        line = line.strip()
        if not line:
            continue

        if line.startswith("&"):
            current_section = line.strip()
            added_keywords[current_section] = {}
        elif line.startswith("/"):
            current_section = None
        elif line.startswith("CELL_PARAMETERS"):
            if line.split()[1] == "angstrom":
                current_section = line.split()[0]
                added_keywords[current_section] = {}
            else:
                print("Only 'CELL_PARAMETERS angstrom' value can be parsed!")
                break
        elif line.startswith("ATOMIC_SPECIES"):
            current_section = line.split()[0]
            added_keywords[current_section] = {}
        elif line.startswith("ATOMIC_POSITIONS"):
            if line.split()[1] == "angstrom":
                current_section = line.split()[0]
                added_keywords[current_section] = {}
            else:
                print("Only 'ATOMIC_POSITIONS angstrom' value can be parsed!")
                break
        elif line.startswith("K_POINTS"):
            if line.split()[1] == "automatic" or line.split()[1] == "gamma":
                current_section = line.split()[0]
                added_keywords[current_section]["K_POINTS"] = line.split()[1]
            else:
                print("Only 'K_POINTS automatic or gamma' value can be parsed!")
                break
        elif current_section:
            n = len(line.split())
            if "=" in line:
                key, value = line.split("=", 1)
                key = key.strip()
                value = value.strip(", ")
                added_keywords[current_section][key] = value
            elif current_section == "ATOMIC_SPECIES" and n == 3:
                species_list = line.split()
                added_keywords[current_section][species_list[0]] = (species_list[1], species_list[2])
            elif current_section == "ATOMIC_POSITIONS" and (n == 7 or n == 4):
                species = line.split()[0]
                coordinates = np.float_(line.split()[1:4])
                if n == 7:
                    constraints = [True if i == '0' else False for i in line.split()[4:7]]
                else:
                    constraints = [False, False, False]
                ase_constraints.append(constraints[0])
                atomic_positions.append(([species, coordinates, constraints]))

            elif current_section == "K_POINTS" and n == 6:
                added_keywords[current_section]["K_POINTS"] += " " + line

    added_keywords["CELL_PARAMETERS"] = atoms.get_cell()
    added_keywords["ATOMIC_POSITIONS"] = atomic_positions

    update_text()
    update_pseudo_buttons()

    atoms.set_constraint(FixAtoms(mask=ase_constraints))
    view(atoms)

def open_qe_out():
    file_path = filedialog.askopenfilename(filetypes=[("QE output files", "*.out")])
    if not file_path:
        return

    atoms = read(file_path, format="espresso-out", index=":")
    view(atoms)

def select_pseudo(element):
    file_path = filedialog.askopenfilename(filetypes=[("Pseudopotential files", "*.UPF")])
    if not file_path:
        return
    file_name = os.path.basename(file_path)
    added_keywords["ATOMIC_SPECIES"][element] = (added_keywords["ATOMIC_SPECIES"][element][0], file_name)
    update_text()

def update_pseudo_buttons():
    for button in pseudo_buttons:
        button.grid_remove()

    row = 6
    for element in added_keywords["ATOMIC_SPECIES"]:
        button = ttk.Button(root, width=15, text=f"Select {element} Pseudo", command=lambda e=element: select_pseudo(e))
        button.grid(column=0, row=row)
        pseudo_buttons.append(button)
        row += 1

def update_keywords(*args):
    keyword_entry["values"] = list(keyword_defaults[section_var.get()].keys())
    keyword_entry.current(0)
    update_input_widget()

def update_text():
    input_text.delete("1.0", tk.END)
    cell_section = True
    for section in sections:
        if section == "&CELL" and not cell_section:
            break
        elif not section.startswith("&"):
            break
        else:
            input_text.insert(tk.END, f"{section}\n")
            for keyword, value in added_keywords[section].items():
                input_text.insert(tk.END, f" {keyword} = {value},\n")
                if keyword == "calculation" and value != "'vc-relax'":
                    cell_section = False
            input_text.insert(tk.END, " /\n")

    if "CELL_PARAMETERS" in added_keywords:
        input_text.insert(tk.END, "CELL_PARAMETERS angstrom\n")
        for row in added_keywords["CELL_PARAMETERS"]:
            input_text.insert(tk.END, f"{row[0]:10.6f} {row[1]:10.6f} {row[2]:10.6f}\n")

    if "ATOMIC_SPECIES" in added_keywords:
        input_text.insert(tk.END, "ATOMIC_SPECIES\n")
        for element, (mass, pseudo_file) in added_keywords["ATOMIC_SPECIES"].items():
            input_text.insert(tk.END, f"{element} {mass} {pseudo_file}\n")

    if "ATOMIC_POSITIONS" in added_keywords:
        input_text.insert(tk.END, "ATOMIC_POSITIONS angstrom\n")
        for element, position, constraint in added_keywords["ATOMIC_POSITIONS"]:
            fixed_status = ['0' if c else '1' for c in constraint]
            input_text.insert(tk.END, f"{element:2s} {position[0]:10.6f} {position[1]:10.6f} {position[2]:10.6f} {' '.join(fixed_status)}\n")

    if "K_POINTS" in added_keywords:
        if added_keywords["K_POINTS"]:
            nk = added_keywords["K_POINTS"]["K_POINTS"].split()
            if nk[0] == "automatic":
                input_text.insert(tk.END, f"K_POINTS {nk[0]}\n")
                input_text.insert(tk.END, f"{nk[1]} {nk[2]} {nk[3]} {nk[4]} {nk[5]} {nk[6]}\n")
            else:
                input_text.insert(tk.END, f"K_POINTS {nk[0]}\n")
        else:
            input_text.insert(tk.END, "K_POINTS automatic\n")
            input_text.insert(tk.END, "1 1 1 0 0 0\n")

def add_keyword():
    section = section_var.get()
    keyword = keyword_var.get()

    if isinstance(keyword_defaults[section][keyword], list):
        value = value_var.get()
    else:
        value = value_entry.get()

    added_keywords[section][keyword] = value
    update_text()

def update_input_widget(*args):
    section = section_var.get()
    keyword = keyword_var.get()
    default_value = keyword_defaults[section][keyword]

    if isinstance(default_value, list):
        value_entry.grid_remove()
        value_dropdown.grid(column=2, row=2)
        value_var.set(default_value[0])
        value_dropdown["values"] = default_value
    else:
        value_dropdown.grid_remove()
        value_entry.grid(column=2, row=2)
        value_entry.delete(0, tk.END)
        value_entry.insert(0, default_value)

def load_cif():
    file_path = filedialog.askopenfilename(filetypes=[("CIF files", "*.cif")])
    if not file_path:
        return

    atoms = read(file_path, format="cif")
    view(atoms)

    added_keywords["CELL_PARAMETERS"] = atoms.get_cell()
    added_keywords["ATOMIC_SPECIES"] = {}
    added_keywords["ATOMIC_POSITIONS"] = {}
    unique_elements = set()

    for atom in atoms:
        element = atom.symbol
        unique_elements.add(element)
        if element not in added_keywords["ATOMIC_SPECIES"]:
            added_keywords["ATOMIC_SPECIES"][element] = (atom.mass, "")

    default_constraints = [[False]*3 for i in range(len(atoms))]
    added_keywords["ATOMIC_POSITIONS"] = [(atom.symbol, atom.position, default_constraints[i]) for i, atom in enumerate(atoms)]
    added_keywords["&SYSTEM"]["ibrav"] = 0
    added_keywords["&SYSTEM"]["nat"] = len(atoms)
    added_keywords["&SYSTEM"]["ntyp"] = len(unique_elements)
    update_text()
    update_pseudo_buttons()

def save_input_file():
    input_text_content = input_text.get("1.0", tk.END)
    file_path = filedialog.asksaveasfilename(defaultextension=".in", filetypes=[("Input files", "*.pwi")])
    if not file_path:
        return
    with open(file_path, "w") as output_file:
        output_file.write(input_text_content)

def show_info():
    messagebox.showinfo("Software Info","QE Input Generator w/ ASE\nby H.-K. Lim\n hklim@kangwon.ac.kr\n")

added_keywords = {section: {} for section in sections}

root = tk.Tk()
root.title("Quantum-Espresso Input Generator")

section_var = tk.StringVar()
section_var.trace("w", update_keywords)
keyword_var = tk.StringVar()
keyword_var.trace("w", update_input_widget)
value_var = tk.StringVar()

load_ase_view_button = ttk.Button(root, text="View QE output", command=open_qe_out, width=15)
load_ase_view_button.grid(column=1, row=0, padx=(0,100))

load_cif_button = ttk.Button(root, text="Load CIF", command=load_cif, width=15)
load_cif_button.grid(column=0, row=1)

section_label = ttk.Label(root, text="Section:")
section_label.grid(column=1, row=0, padx=(200,0))
section_entry = ttk.Combobox(root, width=20, textvariable=section_var)
section_entry["values"] = ["&CONTROL","&SYSTEM","&ELECTRONS","&IONS","&CELL","K_POINTS"]
section_entry.grid(column=2, row=0)

keyword_label = ttk.Label(root, text="Keyword:")
keyword_label.grid(column=1, row=1, padx=(200,0))
keyword_entry = ttk.Combobox(root, width=20, textvariable=keyword_var)
keyword_entry.grid(column=2, row=1)

value_label = ttk.Label(root, text="Value:")
value_label.grid(column=1, row=2, padx=(200,0))
value_entry = ttk.Entry(root, width=20)
value_entry.grid(column=2, row=2)
value_dropdown = ttk.Combobox(root, width=20, textvariable=value_var)
value_dropdown.grid(column=2, row=2)
value_dropdown.grid_remove()

add_button = ttk.Button(root, text="Add Keyword", command=add_keyword, width=15)
add_button.grid(column=2, row=3)

input_text = tk.Text(root, wrap=tk.NONE, width=80, height=30)
input_text.grid(column=0, row=4, columnspan=3)
y_scrollbar = ttk.Scrollbar(root, orient="vertical", command=input_text.yview)
y_scrollbar.grid(column=3, row=4, sticky="ns")
x_scrollbar = ttk.Scrollbar(root, orient="horizontal", command=input_text.xview)
x_scrollbar.grid(column=0, row=5, columnspan=3, sticky="ew")
input_text["yscrollcommand"] = y_scrollbar.set
input_text["xscrollcommand"] = x_scrollbar.set

load_input_button = ttk.Button(root, text="Load Input File", command=load_input_file, width=15)
load_input_button.grid(column=0, row=0)

save_button = ttk.Button(root, text="Save Input File", command=save_input_file, width=15)
save_button.grid(column=0, row=2)

info_button = ttk.Button(root, text="?", command=show_info, width=2)
info_button.grid(column=3,row=6)

pseudo_buttons = []

update_text()

root.mainloop()
