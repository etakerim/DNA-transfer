# Popis:
# Záujmova úloha z biológie - replikácia, translácia a trankripcia
# reťazca DNA. Zahŕňa prepis z DNA na protikus, mRNA, tRNA a priraďi
# triplety k názvom aminokyselín v reťazci. Vstupy normalizuje tak, aby
# bol počet písmen deliteľný tromi (geneticky kód je v tripletoch).
#
# Autor: Miroslav Hájek
# Škola: Gymnázium, Hubeného 23 - sexta (2016/17), septima (2017/18)
# Licencia: GNU GPLv2.1
#
# C-čková verzia programu ------------
# Prezentované: 31.3.2017, Úprava: 27.6.2017
# (Tabuľky aminokyselín a polypeptídové väzby cez SMILES
#  AVOGADRO softvér = zobrazenie 3D modelu)
#
# Python verzia -------
# Prepisovanie dokončené: 27.10.2017
#
# SMILES sú upravené pre ľahké peptídové naviazanie, tj. N je vpredu
# R-COOH + H2N-R, -> R-CO-NH-R,


import sys
import random


triplety = {
    'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAU': 'N',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
    'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGU': 'S',
    'AUA': 'I', 'AUC': 'I', 'AUG': 'M', 'AUU': 'I',
    'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAU': 'H',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
    'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
    'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAU': 'D',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
    'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
    'UAA': ' ', 'UAC': 'Y', 'UAG': ' ', 'UAU': 'Y',
    'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
    'UGA': ' ', 'UGC': 'C', 'UGG': 'W', 'UGU': 'C',
    'UUA': 'L', 'UUC': 'F', 'UUG': 'L', 'UUU': 'F'
}

# Grafické prostredie zvýrazní platné časti dna
aminokyseliny = {
    'A': ('Alanin'              ,'NC(C)C(=O)O'         ),
    ' ': ('STOP'                ,''                    ),
    'C': ('Cystein'             ,'NC(C(=O)O)CS'        ),
    'D': ('Kyselina asparagova' ,'NC(C(=O)O)CC(=O)O'   ),
    'E': ('Kyselina glutamova'  ,'NC(C(=O)O)CCC(=O)O'  ),
    'F': ('Fenylalanin'         ,'NC(C(=O)O)Cc1ccccc1' ),
    'G': ('Glycin'              ,'NCC(=O)O'            ),
    'H': ('Histidin'            ,'NC(C(=O)O)CC1=CN=CN1'),
    'I': ('Izoleucin'           ,'NC(C(=O)O)C(C)CC'    ),
    'K': ('Lyzin'               ,'N(C(=O)O)CCCCN'      ),
    'L': ('Leucin'              ,'NC(C(=O)O)CC(C)C'    ),
    'M': ('Metionin - S'        ,'NC(C(=O)O)CCSC'      ),
    'N': ('Asparagin'           ,'NC(=O)CC(N)C(=O)O'   ),
    'P': ('Prolin'              ,'N1CCCC1C(=O)O'       ),
    'Q': ('Glutamin'            ,'NC(=O)CCC(N)C(=O)O'  ),
    'R': ('Arginin'             ,'NC(CCCNC(N)=N)C(=O)O'),
    'S': ('Serin'               ,'NC(C(=O)O)CO'        ),
    'T': ('Treonin'             ,'NC(C(=O)O)C(O)C'     ),
    'V': ('Valin'               ,'NC(C(=O)O)C(C)C'     ),
    'W': ('Tryptofan'      ,'NC(C(=O)O)CC1=CNc2ccccc12'),
    'Y': ('Tyrozin'           ,'NC(Cc1ccc(O)cc1)C(=O)O')
}


def najdi_doplnok(baza, je_dna):
    baza = baza.upper()

    if baza == 'C':
        return 'G'
    elif baza == 'G':
        return 'C'
    elif baza == 'T' or baza == 'U':
        return 'A'
    elif baza == 'A':
        if je_dna:
            return 'T'
        else:
            return 'U'
    else:
        return ''


def dna_replikacia(zdroj, je_dna):
    ciel = ''
    for b in zdroj:
        ciel += najdi_doplnok(b, je_dna)
    return ciel


def dna_generovat(dna_dlzka):
    nukleotidy = ['A', 'T', 'G', 'C']
    dna_dlzka -= dna_dlzka % 3
    dna = ''

    for i in range(dna_dlzka):
        dna += random.choice(nukleotidy)

    return dna


def dna_skontroluj(dna):
    dna = dna[:len(dna) - (len(dna) % 3)]
    for i, b in enumerate(dna):
        if najdi_doplnok(b, True):
            return i + 1
    return 0


def aminokys_retazec(mRna):
    amino = ''
    for start in range(0, len(mRna), 3):
        kodon = mRna[start:start+3]
        amino += triplety[kodon]

    return amino


def aminokys_nazvy(poly):
    for t in poly:
        yield aminokyseliny[t][0]


# smiles cykly - prečísluj vzory podľa ich získanej pozície
def smiles_cykly(pocitadlo, ak):
    ak_nova = ''
    je_cyklus = False
    for c in aminokyseliny[ak][1]:
        if c == '1' or c == '2':
            t = pocitadlo if c == '1' else pocitadlo + 1
            je_cyklus = True
            if pocitadlo < 10:
                num = '{}'.format(t)
            else:
                num = '%{}'.format(t)
            ak_nova += num
        else:
            ak_nova += c
    if je_cyklus:
        pocitadlo += 1

    return (pocitadlo, ak_nova)


# Skladaj SMILES: Nájdi karboxylovú skupinu C(=O)O
#                 Vymeň O za (<aminokyselina>)
def smiles_vypis(polypeptid):
    inc = 1
    zvysok = ''
    smiles = ''

    for ak in polypeptid:
        inc, ak_nova = smiles_cykly(inc, ak)
        amino = '({}){}'.format(ak_nova, zvysok[1:])

        karboxylO = amino.find("C(=O)O") + len("C(=O)O")
        smiles += amino[:karboxylO - 1]
        zvysok = amino[karboxylO:]

    smiles += zvysok
    return smiles


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
