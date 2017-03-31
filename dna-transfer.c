/*
 * Popis: Zaujmova uloha z biologie - replikacia, translacia a trankripcia
 * retazca DNA. Zahrna prepis z DNA na protikus, mRNA, tRNA a v subore
 * vyhlada nazvy aminokyselin v retazci. Vstupy normalizuje tak, aby 
 * bol pocet pismen delitelny tromi (geneticky kod je v tripletoch).
 *
 * Autor: Miroslav Hajek
 * Škola: Gymnazium, Hubeneho 23 - sexta 2016/2017
 * Licencia: GNU GPLv2.1
 *
 */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>

#define POCET_TRIPL     64
#define TROJICA         3

typedef struct {
    char *d;
    size_t len;
} String;

typedef struct {
    char kodon[5];
    char chem[25];
} Triplet;

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

int tripletcmp(const void *a, const void *b)
{
    return strncmp(((Triplet *)a)->kodon, ((Triplet *)b)->kodon, TROJICA);
}

void aminokys_nacitaj(FILE *fr, Triplet pary[])
{
    String fileline;
    char *token, *token2;
    unsigned int i;

    for (i = 0; i < POCET_TRIPL; i++) {
        fileline = readline(fr);
        if (fileline.d == '\0')
            break;
        
        token = strtok(fileline.d, " ");
        token2 = strtok(NULL, " ");
        if (token != NULL || token2 != NULL) {
            strncpy(pary[i].kodon, token, 5); 
            strncpy(pary[i].chem, token2, 25);
        }
        
        free(fileline.d);
    }
}

size_t aminokys_binhladaj(const char *hl_kodon, Triplet *pary)
{
    int low = 0;
    int high = POCET_TRIPL - 1;
    int dir, mid;

    while (low <= high) {
        mid = (low + high) / 2;
        if ((dir = strncmp(hl_kodon, pary[mid].kodon, TROJICA)) < 0) 
            high = mid - 1;
        else if (dir > 0)
            low = mid + 1;
        else
            return mid;
        }
}

void aminokys_vypis(String mRna, const char *subor, FILE *f)
{
    unsigned int i, findi;
    char kodon[5];
    Triplet triplety[POCET_TRIPL];
    FILE *ak_subor = fopen(subor, "r");
    
    if (ak_subor == NULL) {
        perror("Zoznam aminokyselin nenajdeny");
        user_wait();
        return;
    }
    aminokys_nacitaj(ak_subor, triplety);
    fclose(ak_subor);
    qsort(triplety, POCET_TRIPL, sizeof(Triplet), tripletcmp);

    /* Vylepsene binarnym hladanym z linearneho:
     * O(64n) k O(2(n log 64)) */
    for (i = 0; i < mRna.len; i += TROJICA) {
        memset(kodon, '\0', 5);
        strncpy(kodon, mRna.d + i, TROJICA);
        findi = aminokys_binhladaj(kodon, triplety);

        /* Vypis vysledok */
        if (kodon[0] != '\0' && f != NULL) {
            printf("%s - %s\n", kodon, triplety[findi].chem); 
            fprintf(f, "%s - %s\n", kodon, triplety[findi].chem);
        }
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
    FILE *fw;

    do {
        puts("---- DNA Program ---- ");
        puts("[1]  Nahodna generacia");
        puts("[2]  Zadat vlastny retazec");
        puts("[3]  Zadat cestu suboru");
        printf("Zadajte vas vyber (1 - 3): ");
        vyber = readnumber(stdin);
        putchar('\n');
    } while (vyber < 1 || vyber > 3);

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

    /* Vytvor reťazce DNA2, mRNA, tRNA */
    dna2retazec = dna_replikacia(dna, true);
    m_rna = dna_replikacia(dna, false);
    t_rna = dna_replikacia(m_rna, false);

    
    /* Vypíš výsledky */ 
    fw = fopen("EXPERIMENT.TXT", "w");
    
    fprintf(fw, "DNA1: %s\n", dna.d);
    printf("DNA1: %s\n", dna.d);
    fprintf(fw, "DNA2: %s\n", dna2retazec.d);
    printf("DNA2: %s\n", dna2retazec.d);
    fprintf(fw, "mRNA: %s\n", m_rna.d);
    printf("mRNA: %s\n", m_rna.d);
    fprintf(fw, "tRNA: %s\n", t_rna.d);
    printf("tRNA: %s\n", t_rna.d);
    aminokys_vypis(m_rna, "AMINOKYSELINY.TXT", fw);
    fclose(fw);

    user_wait();

    free(dna.d);
    free(dna2retazec.d);
    free(m_rna.d);
    free(t_rna.d);
} 
