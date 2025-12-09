import os
import numpy as np

def read_fasta(filename):
    """
    Lire un fichier FASTA et retourner une liste de séquences.
    """
    # Obtenir le chemin absolu du fichier
    script_dir = os.path.dirname(os.path.abspath(__file__))
    filepath = os.path.join(script_dir, filename)
    
    sequences = []
    current_seq = ""
    
    with open(filepath, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if current_seq:
                    sequences.append(current_seq)
                    current_seq = ""
            else:
                current_seq += line
        
        if current_seq:
            sequences.append(current_seq)
    
    return sequences


def score_matrix_alignement_global(seq1, seq2, match=1, mismatch=0, gap=-1):
    """
    Créer une matrice de score pour l'alignement global (Needleman-Wunsch).
    Args:
    match: score pour une correspondance (défaut: 1)
    mismatch: score pour une non-correspondance (défaut: 0)
    gap: pénalité pour un gap (défaut: -1)
    """
    n = len(seq1)
    m = len(seq2)
    
    # Créer une matrice de zéros
    matrix = np.zeros((n + 1, m + 1), dtype=int)
    
    # Initialiser la première colonne (gaps dans seq2)
    for i in range(n + 1):
        matrix[i][0] = i * gap
    
    # Initialiser la première ligne (gaps dans seq1)
    for j in range(m + 1):
        matrix[0][j] = j * gap
    
    # Remplir la matrice
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match_score = match if seq1[i-1] == seq2[j-1] else mismatch
            matrix[i][j] = max(
                matrix[i-1][j-1] + match_score,  # Diagonal
                matrix[i-1][j] + gap,            # Haut
                matrix[i][j-1] + gap             # Gauche
            )
    
    return matrix


def score_matrix_alignement_local(seq1, seq2, match=1, mismatch=0, gap=-1):
    """
    Créer une matrice de score pour l'alignement local (Smith-Waterman).
    Returns:
    matrix: la matrice de score
    max_pos: position (i, j) du score maximal pour démarrer le backtracking
    max_score: le score maximal
    """
    n = len(seq1)
    m = len(seq2)
    
    # Créer une matrice de zéros
    matrix = np.zeros((n + 1, m + 1), dtype=int)
    
    # Variables pour suivre le score maximal
    max_score = 0
    max_pos = (0, 0)
    
    # Remplir la matrice
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match_score = match if seq1[i-1] == seq2[j-1] else mismatch
            matrix[i][j] = max(
                0,                                # Minimum = 0 (Spécifique au local)
                matrix[i-1][j-1] + match_score,   # Diagonal
                matrix[i-1][j] + gap,             # Haut
                matrix[i][j-1] + gap              # Gauche
            )
            
            # Mettre à jour le score maximal
            if matrix[i][j] > max_score:
                max_score = matrix[i][j]
                max_pos = (i, j)
    
    return matrix, max_pos, max_score


def backtracking_alignement_global(seq1, seq2, score_matrix, match=1, mismatch=0, gap=-1):
    """
    Backtracking pour l'alignement global.
    Retourne les deux séquences alignées.
    """
    align1 = ""
    align2 = ""
    
    i = len(seq1)
    j = len(seq2)
    
    # Parcourir la matrice de la fin vers le début
    while i > 0 and j > 0:
        score_current = score_matrix[i][j]
        score_diag = score_matrix[i-1][j-1]
        score_up = score_matrix[i-1][j]
        score_left = score_matrix[i][j-1]
        
        match_score = match if seq1[i-1] == seq2[j-1] else mismatch
        
        # Vérifier la provenance du score
        if score_current == score_diag + match_score:
            # Diagonal
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif score_current == score_up + gap:
            # Haut
            align1 = seq1[i-1] + align1
            align2 = "-" + align2
            i -= 1
        elif score_current == score_left + gap:
            # Gauche
            align1 = "-" + align1
            align2 = seq2[j-1] + align2
            j -= 1
        else:
            # Fallback en cas d'erreur d'arrondi ou autre
            i -= 1 
            j -= 1
    
    # Ajouter les gaps restants
    while i > 0:
        align1 = seq1[i-1] + align1
        align2 = "-" + align2
        i -= 1
    
    while j > 0:
        align1 = "-" + align1
        align2 = seq2[j-1] + align2
        j -= 1
    
    return align1, align2


def backtracking_alignement_local(seq1, seq2, score_matrix, max_pos, match=1, mismatch=0, gap=-1):
    """
    Backtracking pour l'alignement local.
    Démarre au score maximal et s'arrête à 0.
    """
    align1 = ""
    align2 = ""
    
    i, j = max_pos  # Commencer à la position du score maximal
    
    # Parcourir jusqu'à atteindre un score de 0
    while i > 0 and j > 0 and score_matrix[i][j] > 0:
        score_current = score_matrix[i][j]
        score_diag = score_matrix[i-1][j-1]
        score_up = score_matrix[i-1][j]
        score_left = score_matrix[i][j-1]
        
        match_score = match if seq1[i-1] == seq2[j-1] else mismatch
        
        if score_current == score_diag + match_score:
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif score_current == score_up + gap:
            align1 = seq1[i-1] + align1
            align2 = "-" + align2
            i -= 1
        elif score_current == score_left + gap:
            align1 = "-" + align1
            align2 = seq2[j-1] + align2
            j -= 1
        else:
            break
    
    return align1, align2


def tableau_alignement(seq1, seq2, matrix):
    """
    Afficher la matrice de manière formatée dans le terminal.
    """
    n = len(seq1)
    m = len(seq2)
    
    # En-tête avec seq2
    print("\n       ", end="")
    print("   -", end="")
    for c in seq2:
        print(f"  {c:>3}", end="")
    print()
    
    # Ligne de séparation
    print("      +" + "-----" * (m + 1))
    
    # Afficher chaque ligne
    for i in range(n + 1):
        if i == 0:
            row_label = "-"
        else:
            row_label = seq1[i - 1]
        
        print(f"   {row_label}  |", end="")
        for j in range(m + 1):
            print(f"  {matrix[i][j]:>3}", end="")
        print()


# --- PROGRAMME PRINCIPAL ---

# Lire les séquences du fichier FASTA
try:
    seqs = read_fasta('sequences.fasta')
    seq1, seq2 = seqs[0], seqs[1]
except Exception as e:
    print(f"Erreur: Assurez-vous que le fichier 'sequences.fasta' existe. ({e})")
    # Valeurs par défaut pour tester si le fichier n'existe pas
    seq1, seq2 = "ATGCGTACGTTAGC", "ATGCCGTCGTTAGG"

print("Séquence 1:", seq1)
print("Séquence 2:", seq2)

# ===================== ALIGNEMENT GLOBAL =====================
print("\n" + "="*60)
print("         ALIGNEMENT GLOBAL (Needleman-Wunsch)")
print("="*60)

# 1. Matrice Global
matrix_global = score_matrix_alignement_global(seq1, seq2)
print("\n=== Matrice de score (Global) ===")
tableau_alignement(seq1, seq2, matrix_global)
print(f"\nScore final: {matrix_global[-1][-1]}")

# 2. Backtracking Global
align1_global, align2_global = backtracking_alignement_global(seq1, seq2, matrix_global)
print("\n=== Alignement Global ===")
print(f"Seq1: {align1_global}")
print(f"Seq2: {align2_global}")

# ===================== ALIGNEMENT LOCAL =====================
print("\n" + "="*60)
print("          ALIGNEMENT LOCAL (Smith-Waterman)")
print("="*60)

# 3. Matrice Local
matrix_local, max_pos, max_score = score_matrix_alignement_local(seq1, seq2)
print("\n=== Matrice de score (Local) ===")
tableau_alignement(seq1, seq2, matrix_local)
print(f"\nScore maximal: {max_score}")
print(f"Position du score maximal: {max_pos}")

# 4. Backtracking Local
align1_local, align2_local = backtracking_alignement_local(seq1, seq2, matrix_local, max_pos)
print("\n=== Alignement Local ===")
print(f"Seq1: {align1_local}")
print(f"Seq2: {align2_local}")