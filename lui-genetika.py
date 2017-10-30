import sys

print('---- DNA Program ----')
print('[1] Náhodné generovanie')
print('[2] Zadat retazec DNA')
print('[3] Zadat DNA zo suboru')
print('[4] Polypeptid skrátene')

vyber = int(input('Zadajte vas vyber (1 - 4): '))

if vyber == 1:
    dna_dlzka = int(input('Zadajte pocet dusikatých báz: '))
    dna = dna_generovat(dna_dlzka)

elif vyber == 2:
    dna = input('Vstup: ')

elif vyber == 3:
    try:
        with open(input('Cela cesta k suboru: '), 'r') as subor:
            dna = subor.read()
    except OSError as err:
        print('Chyba pri otváraní súboru', err)
        sys.exit()

elif vyber == 4:
    polypeptid = input('Vstup Bielkovina: ')
    print(repr(list(aminokys_nazvy(polypeptid))))
    print(smiles_vypis(polypeptid))
    sys.exit()

else:
    print('Neplatný výber')
    sys.exit()

chyba = dna_skontroluj(dna)
if chyba > 1:
    print('Chyba v dusíkatej báze na pozícii {}', chyba)
    sys.exit()

dna_kopia = dna_replikacia(dna, True)
m_rna = dna_replikacia(dna, False)
t_rna = dna_replikacia(m_rna, False)
polypeptid = aminokys_retazec(m_rna)

print('DNA1: {} \n'.format(dna))
print('DNA2: {} \n'.format(dna_kopia))
print('mRNA: {} \n'.format(m_rna))
print('tRNA: {} \n'.format(t_rna))
print('Amino: {} \n'.format(polypeptid))
print(smiles_vypis(polypeptid))
print(repr(list(aminokys_nazvy(polypeptid))))
