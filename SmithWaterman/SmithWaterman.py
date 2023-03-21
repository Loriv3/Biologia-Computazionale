# Author: Lorenzo Tomassetti
#i file gene1.txt e gene2.txt provengono da un parsing dei due geni del file genbank

def SmithWaterman(gene1, gene2):
    # Apertura file contenente il primo gene 
    with open(gene1, "r") as f:
        gene1 = f.read()
        print(gene1)
          
    # Apertura file contenente il secondo gene   
    with open(gene2, "r") as f:
        gene2 = f.read()
        
       
    # Inizializzazione della matrice
    pdMatrix = [[0 for x in range(len(gene1)+1)] for y in range(len(gene2)+1)]
    subMatrix = [[5, -1, -2, -1], [-1, 5, -3, -2], [-2, -3, 5, -2], [-1, -2, -2, 5]]# Matrice di sostituzione
    dictionary = {"A": 0, "C": 1, "G": 2, "T": 3}# Dizionario per accedere alle posizioni della matrice di sostituzione
    maxVal = 0 # Valore massimo della matrice
    maxI = 0 # Posizione i-esima del valore massimo
    maxJ = 0 # Posizione j-esima del valore massimo

    # Riempio matrice secondo regole
    for i in range(1, len(gene2)+1):
        for j in range(1, len(gene1)+1):
            # valori match, delete, insert
            matchOrSub = pdMatrix[i-1][j-1] + subMatrix[dictionary[gene2[i-1]]][dictionary[gene1[j-1]]] 
            delete = pdMatrix[i-1][j] - 4
            insert = pdMatrix[i][j-1] - 4 
            pdMatrix[i][j] = max(0, matchOrSub, delete, insert) # Calcolo valore della cella
            # Controllo se il valore è maggiore del massimo
            # Evita di scorrere di nuovo tutta la matrice
            if maxVal < pdMatrix[i][j]:
                maxVal = pdMatrix[i][j]
                maxI = i
                maxJ = j
            
    
    print("Max value: ", maxVal)
    print("Max value position: ", maxI, maxJ)

    # Traceback
    i = maxI
    j = maxJ
    gene1Aligned = ""
    gene2Aligned = ""
    while pdMatrix[i][j] != 0:
        # Controllo se il valore è stato ottenuto da un match o da una sostituzione
        if pdMatrix[i][j] == pdMatrix[i-1][j-1] + subMatrix[dictionary[gene2[i-1]]][dictionary[gene1[j-1]]]:
            # Match tra le basi
            gene1Aligned = gene1[j-1] + gene1Aligned 
            gene2Aligned = gene2[i-1] + gene2Aligned 
            i -= 1
            j -= 1
        elif pdMatrix[i][j] == pdMatrix[i-1][j] - 4:
            # Gap nel gene1
            gene1Aligned = "-" + gene1Aligned 
            gene2Aligned = gene2[i-1] + gene2Aligned 
            i -= 1
        elif pdMatrix[i][j] == pdMatrix[i][j-1] - 4:
            # Gap nel gene2
            gene1Aligned = gene1[j-1] + gene1Aligned 
            gene2Aligned = "-" + gene2Aligned 
            j -= 1
    print(gene1Aligned) # Stampa gene1 allineato
    print(gene2Aligned) # Stampa gene2 allineato

    

if __name__ == "__main__":
    # Modificare i nomi dei file per fare altri test
    SmithWaterman("gene11.txt", "gene12.txt")