# Genetik
Originaly created as an experimental homework from Biology class with a goal 
of making tidious rewrites of DNA sequence automatic. "Intuitive" command line 
user interface was first developed in C language (still available in legacy/) 
then ported to Python.
This step enabled to create GUI with tkinter -> Genetik.py

### Features
You can generate random DNA sequence of specified length (will be adjusted 
for triplets) or input your own string either on line or from file. There is also 
ability to process protein string to SMILES chemical representation.

###Language
- Slovak (only)

###Screenshot
![alt text](results/screenshot.png?raw=true "In a full glory")


### Command line - Sample execution

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


### Build and Run
Instructions for python version will be posted shortly.
For now: 
```bash
python3 Genetik.py 
```

As simple as compilation can be:
```bash
gcc dna-transfer.c -o program
```
