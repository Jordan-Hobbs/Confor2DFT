import molecules

def main():

    ## gen minimum conformer 
    molecules.find_min_confromer("O=C(OC1=CC(F)=C(C2=CC(F)=C(C#N)C(F)=C2)C(F)=C1)C(C(F)=C3)=CC=C3CCC")

if __name__ == "__main__":
    main()
