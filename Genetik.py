#!/usr/bin/env python3
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from tkinter import font
import gen


class GentikApp(ttk.Frame):

    def __init__(self, root):
        super().__init__(root, padding='2 2 5 5')
        self.grid(sticky='nsew')
        self.root = root
        self.zobraz = []

        root.title('Genetik')
        root.minsize(600, 280)
        root.config(menu=self.menu_creator())
        root.rowconfigure(0, weight=1)
        root.columnconfigure(0, weight=1)

        self.input_fields()
        self.output_configurator()
        self.output_fields()
        self.geometry_manage()

    def menu_creator(self):
        main_menu = tk.Menu(root)

        file_menu = tk.Menu(main_menu, tearoff=0)
        file_menu.add_command(label='Načítaj DNA', command=self.nacitaj_subor)
        file_menu.add_command(label='Ulož pokus', command=self.uloz_pokus)
        file_menu.add_command(label='Ulož SMILES', command=self.uloz_smiles)
        file_menu.add_separator()
        file_menu.add_command(label='Ukonči', command=self.root.quit)

        main_menu.add_cascade(label='Súbor', menu=file_menu)
        return main_menu

    def input_fields(self):
        self.dna = tk.StringVar()
        self.elemcnt = tk.StringVar()
        self.su_tripl = tk.BooleanVar()

        self.in_lbl = ttk.Label(self, text='DNA vstup')
        self.in_entry = ttk.Entry(self, width=50, textvariable=self.dna)
        self.btn_go = ttk.Button(self, text='Prepíš', command=self.sekv_prepis)

        self.pocet_lbl = ttk.Label(self, text='Počet')
        self.num_gen = ttk.Entry(self, width=10, textvariable=self.elemcnt)
        self.bazy_chck = ttk.Radiobutton(self, text='báz',
                                         variable=self.su_tripl, value=False)
        self.trip_chck = ttk.Radiobutton(self, text='tripletov',
                                         variable=self.su_tripl, value=True)
        self.btn_gen = ttk.Button(self, text='Generuj', command=self.sekv_generator)
        self.bazy_chck.invoke()

    def output_configurator(self):
        self.je_dna2 = tk.BooleanVar()
        self.je_mrna = tk.BooleanVar()
        self.je_trna = tk.BooleanVar()
        self.je_ak = tk.BooleanVar()
        self.check_dna2 = ttk.Checkbutton(self, text='DNA2', variable=self.je_dna2,
                                     onvalue=True, offvalue=False)
        self.check_mrna = ttk.Checkbutton(self, text='mRNA', variable=self.je_mrna,
                                     onvalue=True, offvalue=False)
        self.check_trna = ttk.Checkbutton(self, text='tRNA', variable=self.je_trna,
                                     onvalue=True, offvalue=False)
        self.check_ak = ttk.Checkbutton(self, text='AK', variable=self.je_ak,
                                     onvalue=True, offvalue=False)
        self.check_dna2.invoke()
        self.check_mrna.invoke()
        self.check_trna.invoke()
        self.check_ak.invoke()

    def output_fields(self):
        monospace = font.Font(family='Inconsolata', size=12)
        liststyle = {'activestyle': 'none', 'bg': 'white', 'fg': 'black',
                     'bd': 0, 'height': 5, 'font': monospace}

        self.result_lab = tk.Listbox(self, width=6, **liststyle)
        self.result = tk.Listbox(self, width=80, **liststyle)
        self.result.bind('<Double-Button-1>', self.klik_molekula)
        self.resulthand = ttk.Scrollbar(self, orient=tk.HORIZONTAL,
                                        command=self.result.xview)
        self.result.configure(xscrollcommand=self.resulthand.set)
        self.point = ttk.Label(self)

    def geometry_manage(self):
        self.in_lbl.grid(row=0, column=0)
        self.in_entry.grid(row=0, column=1, columnspan=4, pady=10, sticky="we")

        self.check_dna2.grid(row=1, column=1)
        self.check_mrna.grid(row=1, column=2)
        self.check_trna.grid(row=1, column=3)
        self.check_ak.grid(row=1, column=4)
        self.btn_go.grid(row=2, column=2, columnspan=2, sticky="we", pady=10)

        self.pocet_lbl.grid(row=0, column=6, sticky="w")
        self.num_gen.grid(row=0, column=7, sticky="we")
        self.bazy_chck.grid(row=1, column=6, sticky="we")
        self.trip_chck.grid(row=1, column=7, sticky="we")
        self.btn_gen.grid(row=2, column=6,  columnspan=2, sticky="we")

        self.result_lab.grid(row=3, column=0)
        self.result.grid(row=3, column=1, columnspan=8)
        self.resulthand.grid(row=4, column=1, columnspan=8, sticky="we")
        self.point.grid(row=7, column=2, columnspan=3, sticky="w")

    def sekv_generator(self):
        try:
            n = int(self.elemcnt.get())
            if self.su_tripl.get():
                n *= 3
            self.dna.set(gen.dna_generovat(n))
        except:
            pass

    def nacitaj_subor(self):
        nazov = filedialog.askopenfilename()
        if nazov:
            with open(nazov, 'r') as subor:
                self.dna.set(subor.read())

    def uloz_pokus(self):
        subor = filedialog.asksaveasfile(mode='w', defaultextension='.txt')
        if subor:
            for polozka in self.zobraz:
                subor.write('{}, {}\n'.format(polozka[0], polozka[1]))
            subor.close()

    def uloz_smiles(self):
        subor = filedialog.asksaveasfile(mode='w', defaultextension='.smi')
        if subor:
            ak = gen.aminokys_retazec(gen.dna_replikacia(self.dna.get()))
            gen.smiles_vypis(subor, ak)
            subor.close()

    def klik_molekula(self, click):
        i_active = self.result.index(tk.ACTIVE)
        bx, by, bw, bh = self.result.bbox(i_active)
        s = self.zobraz[i_active][1]
        typ = self.zobraz[i_active][0]

        index = (click.x + abs(bx)) // (bw // len(s))
        if index < len(s):
            if typ == 'Amino':
                popis = gen.aminokyseliny[s[index]][0]
            else:
                popis = gen.bazy[s[index]]
        self.point['text'] = 'Látka: {} ({})'.format(popis, s[index])

    def ak_padding(self, ak):
        ak_pad = ''
        for p in ak:
            if p.isspace():
                ak_pad += ' ' * 3
            else:
                ak_pad += '-{}-'.format(p)
        return ak_pad

    def sekv_prepis(self):
        self.result_lab.delete(0, tk.END)
        self.result.delete(0, tk.END)

        dna = self.dna.get()
        m_rna = gen.dna_replikacia(dna)

        self.zobraz = [('DNA-A', dna)]
        if self.je_dna2.get():
            self.zobraz.append(('DNA-B', gen.dna_replikacia(dna, True)))
        if self.je_mrna.get():
            self.zobraz.append(('mRNA', m_rna))
        if self.je_trna.get():
            self.zobraz.append(('tRNA', gen.dna_replikacia(m_rna)))
        if self.je_ak.get():
            self.zobraz.append(('Amino', gen.aminokys_retazec(m_rna)))

        for item in self.zobraz:
            self.result_lab.insert(tk.END, item[0])
            if item[0] == 'Amino':
                self.result.insert(tk.END, self.ak_padding(item[1]))
            else:
                self.result.insert(tk.END, item[1])


if __name__ == '__main__':
    root = tk.Tk()
    app = GentikApp(root)
    app.mainloop()
