from tkinter import *
from tkinter import filedialog
from tkinter import ttk

def dna_generuj():
    import random
    try:
        n = int(elemcnt.get())
        if not su_bazy.get():
            n *= 3
        s = ''
        for i in range(n):
            s += random.choice(['A', 'T', 'G', 'C'])
        dna.set(s)
    except:
        pass

def nacitaj_subor():
    nazov = filedialog.askopenfilename()
    if not nazov:
        return
    with open(nazov, 'r') as subor:
        dna.set(subor.read())

def uloz_pokus():
    subor = filedialog.asksaveasfile(mode='w', defaultextension='.txt')
    if subor:
        subor.write(dna.get())
        subor.close()

result = None
point = None

def elem_active(click):
    active = result.index(ACTIVE)
    bx, by, bw, bh = result.bbox(active)
    s = result.get(active)

    index = click.x // (bw // len(s))
    if index < len(s):
        point['text'] = s[index]

root = Tk()
root.title('Genetika')
root.minsize(600, 250)


main_menu = Menu(root)

file_menu = Menu(main_menu, tearoff=0)
file_menu.add_command(label='Načítaj DNA', command=nacitaj_subor)
file_menu.add_command(label='Ulož pokus', command=uloz_pokus)
file_menu.add_command(label='Ulož SMILES', command=uloz_pokus)
file_menu.add_separator()
file_menu.add_command(label='Ukonči', command=root.quit)

main_menu.add_cascade(label='Súbor', menu=file_menu)
root.config(menu=main_menu)


frame = ttk.Frame(root, padding='2 2 5 5')

dna = StringVar()
elemcnt = StringVar()
su_bazy = BooleanVar()

lab_dna = ttk.Label(frame, text="DNA")
entry_dna = ttk.Entry(frame, width=50, textvariable=dna)
btn_go = ttk.Button(frame, text="Prepíš")

count_gen = ttk.Entry(frame, width=10, textvariable=elemcnt)
bazy_box = ttk.Radiobutton(frame, text="báz",
                           variable=su_bazy, value=True)
trip_box = ttk.Radiobutton(frame, text="tripletov",
                           variable=su_bazy, value=False)
btn_gen = ttk.Button(frame, text="Generuj", command=dna_generuj)
bazy_box.invoke()

je_dna2 = BooleanVar()
je_mrna = BooleanVar()
je_trna = BooleanVar()
je_ak = BooleanVar()
check_dna2 = ttk.Checkbutton(frame, text="DNA2", variable=je_dna2,
                             onvalue=True, offvalue=False)
check_mrna = ttk.Checkbutton(frame, text="mRNA", variable=je_mrna,
                             onvalue=True, offvalue=False)
check_trna = ttk.Checkbutton(frame, text="tRNA", variable=je_trna,
                             onvalue=True, offvalue=False)
check_ak = ttk.Checkbutton(frame, text="AK", variable=je_ak,
                             onvalue=True, offvalue=False)
check_dna2.invoke()
check_mrna.invoke()
check_trna.invoke()
check_ak.invoke()

result_lab = Listbox(frame, width=5, height=5)
result_lab.insert(END, "DNA1")

result = Listbox(frame, width=80, height=5)
result.bind('<Button-1>', elem_active)
result.insert(END, 'ATCGGGATCCCTGACAGATCAGTACGTTTGACGAAATGACCCAGTATTG')
result.insert(END, 'AHOJ')

resulthand = ttk.Scrollbar(frame, orient=HORIZONTAL, command=result.xview)
result.configure(xscrollcommand=resulthand.set)

point = ttk.Label(frame)

frame.grid(row=0, column=0, sticky="nsew")

lab_dna.grid(row=0, column=0)
entry_dna.grid(row=0, column=1, columnspan=5, pady=10)
check_dna2.grid(row=1, column=1)
check_mrna.grid(row=1, column=2)
check_trna.grid(row=1, column=3)
check_ak.grid(row=1, column=4)
btn_go.grid(row=2, column=2, columnspan=3)

btn_gen.grid(row=2, column=7,  columnspan=8)
count_gen.grid(row=0, column=7, columnspan=8)
bazy_box.grid(row=1, column=7)
trip_box.grid(row=1, column=8)

result_lab.grid(row=4, column=0)
result.grid(row=4, column=1, columnspan=8)
resulthand.grid(row=5, column=1, columnspan=8, sticky="we")
point.grid(row=7, column=2)

#frame.rowconfigure(0, weight=1)
root.rowconfigure(0, weight=1)
root.columnconfigure(0, weight=1)

root.mainloop()
