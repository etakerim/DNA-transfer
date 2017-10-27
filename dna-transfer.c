/*
 * Popis: Zaujmova uloha z biologie - replikacia, translacia a trankripcia
 * retazca DNA. Zahrna prepis z DNA na protikus, mRNA, tRNA a v subore
 * vyhlada nazvy aminokyselin v retazci. Vstupy normalizuje tak, aby
 * bol pocet pismen delitelny tromi (geneticky kod je v tripletoch).
 *
 * Autor: Miroslav Hajek
 * Škola: Gymnazium, Hubeneho 23 - sexta 2016/2017
 * Licencia: GNU GPLv2.1
 * Prezentované: 31.3.2017
 * Úprava:       27.6.2017
 * (Tabuľky aminokyselín a polypeptídové väzby cez SMILES
 *  AVOGADRO softvér = zobrazenie 3D modelu)
 */

/* SMILES sú upravené pre ľahké peptídové naviazanie, tj. N je vpredu
 * R-COOH + H2N-R, -> R-CO-NH-R, */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>

#define POCET_TRIPL     64
#define TROJICA         3
#define STOP            ' '

static const struct {
    char id;
    char *chemnazov;
    char *smiles;
} aminokyseliny[] = {
    {'A', "Alanin"              ,"NC(C)C(=O)O"         },
    {'B', "STOP"                ,""                    },
    {'C', "Cystein"             ,"NC(C(=O)O)CS"        },
    {'D', "Kyselina asparagova" ,"NC(C(=O)O)CC(=O)O"   },
    {'E', "Kyselina glutamova"  ,"NC(C(=O)O)CCC(=O)O"  },
    {'F', "Fenylalanin"         ,"NC(C(=O)O)Cc1ccccc1" },
    {'G', "Glycin"              ,"NCC(=O)O"            },
    {'H', "Histidin"            ,"NC(C(=O)O)CC1=CN=CN1"},
    {'I', "Izoleucin"           ,"NC(C(=O)O)C(C)CC"    },
    {'J', ""                    ,""                    },
    {'K', "Lyzin"               ,"N(C(=O)O)CCCCN"      },
    {'L', "Leucin"              ,"NC(C(=O)O)CC(C)C"    },
    {'M', "Metionin"            ,"NC(C(=O)O)CCSC"      },
    {'N', "Asparagin"           ,"NC(=O)CC(N)C(=O)O"   },
    {'O', ""                    ,""                    },
    {'P', "Prolin"              ,"N1CCCC1C(=O)O"       },
    {'Q', "Glutamin"            ,"NC(=O)CCC(N)C(=O)O"  },
    {'R', "Arginin"             ,"NC(CCCNC(N)=N)C(=O)O"},
    {'S', "Serin"               ,"NC(C(=O)O)CO"        },
    {'T', "Treonin"             ,"NC(C(=O)O)C(O)C"     },
    {'U', ""                    ,""                    },
    {'V', "Valin"               ,"NC(C(=O)O)C(C)C"     },
    {'W', "Tryptofan"      ,"NC(C(=O)O)CC1=CNc2ccccc12"},
    {'X', ""                    ,""                    },
    {'Y', "Tyrozin"           ,"NC(Cc1ccc(O)cc1)C(=O)O"}
};

static const struct {
    char *kodon;
    char id;
} triplety[POCET_TRIPL] = {
    {"AAA", 'K'}, {"AAC", 'N'}, {"AAG", 'K'}, {"AAU", 'N'},
    {"ACA", 'T'}, {"ACC", 'T'}, {"ACG", 'T'}, {"ACU", 'T'},
    {"AGA", 'R'}, {"AGC", 'S'}, {"AGG", 'R'}, {"AGU", 'S'},
    {"AUA", 'I'}, {"AUC", 'I'}, {"AUG", 'M'}, {"AUU", 'I'},
    {"CAA", 'Q'}, {"CAC", 'H'}, {"CAG", 'Q'}, {"CAU", 'H'},
    {"CCA", 'P'}, {"CCC", 'P'}, {"CCG", 'P'}, {"CCU", 'P'},
    {"CGA", 'R'}, {"CGC", 'R'}, {"CGG", 'R'}, {"CGU", 'R'},
    {"CUA", 'L'}, {"CUC", 'L'}, {"CUG", 'L'}, {"CUU", 'L'},
    {"GAA", 'E'}, {"GAC", 'D'}, {"GAG", 'E'}, {"GAU", 'D'},
    {"GCA", 'A'}, {"GCC", 'A'}, {"GCG", 'A'}, {"GCU", 'A'},
    {"GGA", 'G'}, {"GGC", 'G'}, {"GGG", 'G'}, {"GGU", 'G'},
    {"GUA", 'V'}, {"GUC", 'V'}, {"GUG", 'V'}, {"GUU", 'V'},
    {"UAA", ' '}, {"UAC", 'Y'}, {"UAG", ' '}, {"UAU", 'Y'},
    {"UCA", 'S'}, {"UCC", 'S'}, {"UCG", 'S'}, {"UCU", 'S'},
    {"UGA", ' '}, {"UGC", 'C'}, {"UGG", 'W'}, {"UGU", 'C'},
    {"UUA", 'L'}, {"UUC", 'F'}, {"UUG", 'L'}, {"UUU", 'F'}
};

typedef struct {
    char *d;
    size_t len;
} String;

void user_wait(void)
{
    printf("\nStlacte ENTER na ukoncenie...\n");
    getchar();
}

String string_constalloc(size_t len)
{
    String str;

    str.d = malloc(sizeof(char) * (len + 1));
    if (str.d == NULL) {
        fprintf(stderr, "Chyba - Nedostatok pamate RAM\n");
        user_wait();
        exit(3);
    }
    str.len = len;
    str.d[str.len] = '\0';

    return str;
}

String readline(FILE *fr)
{
    char z;
    int allocsz = 8;
    int index = 0;
    String line = string_constalloc(allocsz);

    while ((z = fgetc(fr)) != EOF && z != '\n') {
        line.d[index++] = z;
        if (index == allocsz) {
            allocsz *= 1.5;
            line.d = realloc(line.d, allocsz);
        }
    }

    line.d[index] = '\0';
    line.len = index;
    return line;
}

int readnumber(FILE *fr)
{
    String vstup = readline(fr);
    int cislo = atoi(vstup.d);

    free(vstup.d);
    return cislo;
}

String readtxtfile(void)
{
    String filestr;
    String cesta = readline(stdin);
    FILE *subor = fopen(cesta.d, "r");

    if (subor == NULL) {
        free(cesta.d);
        perror("Chyba pri otvarani suboru");
        user_wait();
        exit(2);
    }

    filestr = readline(subor);
    free(cesta.d);
    fclose(subor);

    return filestr;
}

FILE *subor_zapis_volny(char *prefix, char *sufix)
{
    char filename[256];
    int cnt = 0;
    FILE *fr = NULL;

    do {
        if (fr != NULL)
            fclose(fr);
        snprintf(filename, sizeof(filename), "%s%d.%s", prefix, cnt++, sufix);
        fr = fopen(filename, "r");
    } while (fr != NULL);

    return fopen(filename, "w");
}

char najdi_doplnok(char baza, bool je_dna)
{
    baza = toupper(baza);

    switch (baza) {
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        case 'U': return 'A';
        case 'A':
            if (je_dna)
                return 'T';
            else
                return 'U';
        default:
            return '\0';
    }
}

void dna_skrat_podieltri(String *retazec)
{
    while (retazec->len % TROJICA != 0)
        retazec->len--;
    retazec->d[retazec->len] = '\0';
}

int dna_skontroluj(String retazec)
{
    unsigned int i;

    dna_skrat_podieltri(&retazec);
    for (i = 0; i < retazec.len; i++)
        if (najdi_doplnok(retazec.d[i], true) == '\0')
            return i + 1;

    return 0;
}

String dna_generovat(unsigned int dnadlzka)
{
    unsigned int i;
    int rndpos;
    static const char nukleotidy[] = {'A', 'T', 'G', 'C'};
    String dna_retazec = string_constalloc(dnadlzka);

    dna_skrat_podieltri(&dna_retazec);
    srand(time(NULL));
    for (i = 0; i < dna_retazec.len; i++) {
        rndpos = rand() % 4;
        dna_retazec.d[i] = nukleotidy[rndpos];
    }

    return dna_retazec;
}

String dna_replikacia(String zdroj, bool je_dna)
{
    unsigned int i;
    String ciel = string_constalloc(zdroj.len);

    for (i = 0; i < ciel.len; i++) {
        ciel.d[i] = najdi_doplnok(zdroj.d[i], je_dna);
    }

    return ciel;
}

unsigned int hash_ohodnot(char znak, unsigned int index)
{
    const unsigned int faktor[] = {16, 4, 1};
    unsigned int cena;

    switch (znak) {
        case 'A':
            cena = 0; break;
        case 'C':
            cena = 1; break;
        case 'G':
            cena = 2; break;
        case 'U':
            cena = 3; break;
    }

    return cena * faktor[index];
}

int idamino_to_index(char id)
{
    id = toupper(id);
    if ((id < 'A' || id > 'Y') && id != ' ') {
        printf("\nChyba - Neplatne oznacenie aminokyseliny: %c\n", id);
        user_wait();
        exit(2);
    }
    return (id != STOP) ? id - 'A':  1;
}

String aminokys_retazec(String mRna)
{
    unsigned int i, j, index;
    unsigned int ia = 0;
    char kodon[6];
    String amino = string_constalloc(mRna.len / 3);

    for (i = 0; i < mRna.len; i += TROJICA) {
        index = 0;
        memset(kodon, 0, sizeof(kodon));
        strncpy(kodon, mRna.d + i, TROJICA);
        if (kodon[0] == '\0')
            break;

        // Nájdi kodón v abecedne zoradenej tabuľke
        for (j = 0; j < TROJICA; j++)
            index += hash_ohodnot(kodon[j], j);

        // Pridaj symbol do výsledného reťazca
        amino.d[ia++] = triplety[index].id;
    }

    return amino;
}

void aminokys_nazvy(FILE *fw, String poly)
{
    int i, id;

    for (i = 0; i < poly.len; i++) {
        id = idamino_to_index(poly.d[i]);
        fprintf(fw, "%s, ", aminokyseliny[id].chemnazov);
    }
    fprintf(fw, "\n\n");
}


// Pri cykloch v SMILES inkrementovať počítanie
static const char *SMILES_cykly(char *amino)
{
    static int inc = 1;
    static char curr[64];
    char tmp[64];
    char num[32];
    bool jecyklus = false;
    int i, j, k, t;
    int len = strlen(amino);

    strncpy(tmp, amino, sizeof(tmp));
    for (i = 0, j = 0; i < len || j < sizeof(tmp); i++) {
        if (tmp[i] == '1' || tmp[i] == '2') {
            t = (tmp[i] == '1') ? inc : inc + 1;
            jecyklus = true;

            if (inc < 10)
                snprintf(num, sizeof(num), "%d", t);
            else
                snprintf(num, sizeof(num), "%%%d", t);
            for (k = 0; num[k] != '\0'; k++)
                curr[j++] = num[k];

        } else {
            curr[j++] = tmp[i];
        }
    }

    if (jecyklus)
        inc++;
    return curr;
}

void SMILES_polypeptid(FILE *dst, String polypeptid)
{
    /* Skladaj SMILES: Nájdi karboxylovú skupinu C(=O)O
                       Vymeň O za (<aminokyselina>)
        zvysokprev = ""
        cyklus:
             amino  = '(' + amino ')' + zvysokprev
             i = smiles.find("C(=O)O") + len("C(=O)O")
             save = smiles[:i-1]   //print
             zvysokprev = smiles[i:]
         print zvysok
    */
    #define BUF_LEN     1024
    int i, pos;
    char amino[BUF_LEN] = {};
    char zvysok[BUF_LEN] = {};
    char ulozit[BUF_LEN] = {};
    char *karboxylO = NULL;

    for (i = 0; i < polypeptid.len; i++) {
        pos = idamino_to_index(polypeptid.d[i]);

        memset(amino, 0, sizeof(amino));
        snprintf(amino, sizeof(amino),
                "(%s)%s", SMILES_cykly(aminokyseliny[pos].smiles), zvysok + 1);

        karboxylO = strstr(amino, "C(=O)O");
        if (karboxylO != NULL) {
            karboxylO += 5;
            memset(ulozit, 0, sizeof(ulozit));
            strncpy(ulozit, amino, karboxylO - amino);
            fprintf(dst, "%s", ulozit);

            memset(zvysok, 0, sizeof(zvysok));
            strcpy(zvysok, karboxylO);
        }
    }

    fprintf(dst, "%s\n", zvysok);
}

void SMILES_subor(String polypeptid, char *prefix)
{
    FILE *fsmile;

    if (polypeptid.len > 333) {
        fprintf(stderr, ("SMILES Nebude generovany!: pocet"
                        "prvkov proteinu je %d z max. povolenych 333\n"), polypeptid.len);
        return;
    }

    fsmile = subor_zapis_volny(prefix, "smiles");
    if (fsmile != NULL) {
        SMILES_polypeptid(fsmile, polypeptid);
        fclose(fsmile);
    }
}


int main(void)
{
    int vyber;
    unsigned int dlzkadna;
    String dna;
    String dna2retazec;
    String m_rna;
    String t_rna;
    String polypeptid;
    FILE *fw;

    do {
        puts("---- DNA Program ---- ");
        puts("[1]  Nahodna generacia");
        puts("[2]  Zadat vlastny retazec");
        puts("[3]  Zadat cestu suboru");
        puts("[4]  Polypeptid skratene");
        printf("Zadajte vas vyber (1 - 4): ");
        vyber = readnumber(stdin);
        putchar('\n');
    } while (vyber < 1 || vyber > 4);

    switch (vyber) {
        case 1:
            printf("Zadajte pocet dusikatych baz: ");
            dlzkadna = readnumber(stdin);
            dna = dna_generovat(dlzkadna);
            break;
        case 2:
            printf("Vstup: ");
            dna = readline(stdin);
            break;
        case 3:
            printf("Cela cesta k suboru: ");
            dna = readtxtfile();
            break;
        case 4:
            printf("Vstup Bielkovina: ");
            polypeptid = readline(stdin);
            aminokys_nazvy(stdout, polypeptid);
            SMILES_subor(polypeptid, "bielkovina");
            user_wait();
            exit(0);
        default:
            puts("Neplatny vyber");
            user_wait();
            exit(1);
    }

    if ((vyber = dna_skontroluj(dna)) != 0) {
        printf("Chyba v dusikatej baze DNA na pozicii: %d\n", vyber);
        user_wait();
        exit(3);
    }

    // Vytvor reťazce DNA2, mRNA, tRNA, aminokyseliny
    dna2retazec = dna_replikacia(dna, true);
    m_rna = dna_replikacia(dna, false);
    t_rna = dna_replikacia(m_rna, false);
    polypeptid = aminokys_retazec(m_rna);

    printf("DNA1: %s\n\n", dna.d);
    printf("DNA2: %s\n\n", dna2retazec.d);
    printf("mRNA: %s\n\n", m_rna.d);
    printf("tRNA: %s\n\n", t_rna.d);
    printf("Amino: %s\n\n", polypeptid.d);


    // Vypíš výsledky
    fw = subor_zapis_volny("Pokus", "txt");
    if (fw != NULL) {
        fprintf(fw, "DNA1: %s\n\n", dna.d);
        fprintf(fw, "DNA2: %s\n\n", dna2retazec.d);
        fprintf(fw, "mRNA: %s\n\n", m_rna.d);
        fprintf(fw, "tRNA: %s\n\n", t_rna.d);
        fprintf(fw, "Amino: %s\n\n", polypeptid.d);
        aminokys_nazvy(fw, polypeptid);
        fclose(fw);
    }

    SMILES_subor(polypeptid, "protein");
    free(polypeptid.d);
    free(dna.d);
    free(dna2retazec.d);
    free(m_rna.d);
    free(t_rna.d);

    user_wait();
}
