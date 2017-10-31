# Genetik
Python Tkinter application for demostration of DNA replication, translation
and transcription to RNA and then to Amino acid chain.

### About  
Originaly created as an experimental homework from Biology class with a goal 
of making tidious rewrites of DNA sequence automatic. "Intuitive" command line 
user interface was first developed in C language (still available in legacy/) 
then ported to Python. This step enabled to create GUI with tkinter -> Genetik.py

### Screenshot
![alt text](results/screenshot.png?raw=true "Genetik app")


### Features
- Input your own DNA sequence (A, T, C, G) from within app or from file
- Generate random DNA of specified length
- Export whole "experiment" to file and amino acids to .smiles
- Basic information about meaning of characters


### Installation
Prerequisites: git, python3, pip3

##### Developer build for Users
To distibute it to others consider using `pyinstaller`.
```bash
$ cd Genetik/Genetik-app
$ pip install pyinstaller
$ pyinstaller -w -i genetik.ico -F Genetik.py
# pyinstaller --noconsole --onefile -i <icon=genetik.ico> <app=Genetik.py>
```

##### Development
```bash
$ git clone https://github.com/etakerim/Genetik
$ cd Genetik/Genetik-app
$ python Genetik.py 
```
If you wish to run old command line version
```bash
$ cd Genetik/legacy
$ gcc dna-transfer.c -o genetik
$ ./genetik
```

### How to help?
Basic motivation was to create something that could demonstrate use
of IT in Biology, especially in Genetics. After this first draft
it would be nice if we could make app more captive for students.
If you feel like you want to help, here are some things to get you
started:
- Translations (English, German, ...)
- Better graphical style
- Grafical info about chemicals in DNA and in Amino Acids
- Mendel laws of inheritance
