"""
Modul gen.py obsahuje tabuľky priradení tripletov na aminokyselinu
a aminokyselín na ich názvy. Funkcie na generovanie, replikáciu,
transláciu a transkripciu DNA na všetky ostatné formy DNA2, mRNA, tRNA.
Podpora formátu .dna (čistá sekvencia) a vytváranie .smiles (.smi)


Vytvorené ako záujmova úloha z biológie.
Autor: Miroslav Hájek
Škola: Gymnázium, Hubeného 23
Trieda: sexta, septima (2016-18)
Licencia: GNU LGPLv2.1
"""
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
    'M': ('Metionin'            ,'NC(C(=O)O)CCSC'      ),
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

bazy = {
    'A': 'Adenín',
    'G': 'Guanín',
    'C': 'Cytozín',
    'U': 'Uracil',
    'T': 'Tymín'
}


def najdi_doplnok(baza, repl=False):
    ''' Vráti komplementárny pár k platnej dusíkatej báze 'baza'.
    Ak sa jedná o replikáciu je potrebné zadat vlajku 'repl'.'''
    baza = baza.upper()

    if baza == 'C':
        return 'G'
    elif baza == 'G':
        return 'C'
    elif baza == 'T' or baza == 'U':
        return 'A'
    elif baza == 'A':
        if repl:
            return 'T'
        else:
            return 'U'
    else:
        return ''


def dna_replikacia(zdroj, repl=False):
    '''Vytvorí prepis pre celý retazec dna/rna 'zdroj' s využitím
    'najdi_doplnok' '''
    ciel = ''
    for b in zdroj:
        ciel += najdi_doplnok(b, repl)
    return ciel


def dna_generovat(nbaz):
    ''' Generuje DNA sekvenciu pozostávajúcu zo zadaného počtu
    dusíkatých báz 'nbaz'. Ak nie je počet tripletový zabezpečí to.'''
    nukleotidy = ['A', 'T', 'G', 'C']
    nbaz -= nbaz % 3
    dna = ''

    for i in range(nbaz):
        dna += random.choice(nukleotidy)

    return dna


def dna_skontroluj(dna):
    '''Skontroluje dna sekvenciu. Pokiaľ je nesprávna
    vráti neplatné miesto počnúc od 1. 0 znamená ok sekvenciu '''
    if len(dna) % 3 != 0:
       return 1

    for i, b in enumerate(dna):
        if not najdi_doplnok(b, True):
            return i + 1
    return 0


def aminokys_retazec(mRna):
    '''Z platného reťazca mRNA vytvorí reťazec značiek aminokyselín'''
    amino = ''
    for start in range(0, len(mRna), 3):
        kodon = mRna[start:start+3]
        if len(kodon) == 3:
            amino += triplety[kodon]

    return amino


def aminokys_nazvy(poly):
    for t in poly:
        yield aminokyseliny[t][0]


def smiles_cykly(pocitadlo, ak):
    '''Prečísluje cykly v tabuľkových SMILES podľa
    ich získanej pozície v reťazci'''
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


def smiles_vypis(polypeptid):
    '''Skladá SMILES všetkých AK z polypeptidového reťazca dokopy.
    Nájde karboxylovú skupinu C(=O)Oa vymení  O za (<aminokyselina>) '''
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
