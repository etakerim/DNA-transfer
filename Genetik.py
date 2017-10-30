import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
import gen


class GentikApp(ttk.Frame):

    def __init__(self, root):
        super().__init__(root, padding='2 2 5 5')
        super().grid(sticky='nsew')
        self.root = root
        root.title('Genetika')
        root.minsize(600, 250)
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
        file_menu.add_command(label='Ulož SMILES', command=self.uloz_pokus)
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
        self.btn_go = ttk.Button(self, text='Prepíš')

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
        self.result_lab = tk.Listbox(self, width=5, height=5)
        self.result_lab.insert(tk.END, "DNA1")

        self.result = tk.Listbox(self, width=80, height=5)
        self.result.bind('<Button-1>', self.elem_active)
        self.result.insert(tk.END, 'ATCGGGATCCCTGACAGATCAGTACGTTTGACGAAATGACCCAGTATTG')
        self.result.insert(tk.END, 'AHOJ')

        self.resulthand = ttk.Scrollbar(self, orient=tk.HORIZONTAL,
                                        command=self.result.xview)
        self.result.configure(xscrollcommand=self.resulthand.set)
        self.point = ttk.Label(self)

    def geometry_manage(self):
        self.in_lbl.grid(row=0, column=0)
        self.in_entry.grid(row=0, column=1, columnspan=5, pady=10)
        self.check_dna2.grid(row=1, column=1)
        self.check_mrna.grid(row=1, column=2)
        self.check_trna.grid(row=1, column=3)
        self.check_ak.grid(row=1, column=4)
        self.btn_go.grid(row=2, column=2, columnspan=3)

        self.btn_gen.grid(row=2, column=7,  columnspan=8)
        self.num_gen.grid(row=0, column=7, columnspan=8)
        self.bazy_chck.grid(row=1, column=7)
        self.trip_chck.grid(row=1, column=8)

        self.result_lab.grid(row=4, column=0)
        self.result.grid(row=4, column=1, columnspan=8)
        self.resulthand.grid(row=5, column=1, columnspan=8, sticky="we")
        self.point.grid(row=7, column=2)

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
            subor.write(self.dna.get())
            subor.close()

    def elem_active(self, click):
        active = self.result.index(tk.ACTIVE)
        bx, by, bw, bh = self.result.bbox(active)
        s = self.result.get(active)

        index = click.x // (bw // len(s))
        if index < len(s):
            self.point['text'] = s[index]


if __name__ == '__main__':
    root = tk.Tk()
    app = GentikApp(root)
    app.mainloop()
