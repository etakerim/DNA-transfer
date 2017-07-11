# DNA transfer
Biology class experimental homework of making tidious rewrites and lookups automatic.
Through "intuitive" command line user interface you can generate random DNA sequence
of specified length (adjusted for triplets) or input your own string either
on line or from file. There is also ability to process protein string to
SMILES chemical representation:


### Sample execution

```bash
---- DNA Program ---- 
[1]  Nahodna generacia
[2]  Zadat vlastny retazec
[3]  Zadat cestu suboru
[4]  Polypeptid skratene
Zadajte vas vyber (1 - 4): 1

Zadajte pocet dusikatych baz: 9
DNA1: CTCCATCTC

DNA2: GAGGTAGAG

mRNA: GAGGUAGAG

tRNA: CUCCAUCUC

Amino: EVE


Stlacte ENTER na ukoncenie...
```


### Generated files
After each "experiment" two files will be generated in place of executable:
- Pokus.txt         -> Save generated data from command line 
- protein.smiles    -> SMILES string of protein (Maybe you have to strip leading
                       and trailing paranthesis, but it is valid)

### Build and Run
As simple as compilation can be:
```bash
gcc dna-transfer.c -o program
```

