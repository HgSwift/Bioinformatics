import numpy as np

def optimal_global_alignment(path_fasta_1, path_fasta_2, path_score_mat, gap_penalty = -200):
    score = None
    aligned_1 = ""
    aligned_2 = ""
    f1 = open(path_fasta_1, "r")
    seq1 =f1.readline(6)
    seq1 =f1.read()
    sequence1 = seq1[0:(len(seq1)-1)]
    f2 = open(path_fasta_2, "r")
    seq2 =f2.readline(6)
    seq2 =f2.read()
    sequence2 = seq2[0:(len(seq2)-1)]
    #Read in matrix
    m = open(path_score_mat, "r")
    mat = m.readline()
    mat = m.read()
    sub_m  = [[]]
    sub_row = []
    sub_cell = ""
    for c in mat:
        if c != 'A' and c != 'G' and c != 'T' and c != 'C' and c != '\n':
            if c == '\t':
                sub_row.append(sub_cell)
                sub_cell = ""
            else:
                sub_cell +=c
        elif c == '\n':
            sub_row.append(sub_cell)
            sub_m.append(sub_row[1:])
            sub_row = []
            sub_cell = ""

    sub_m = sub_m[1:]
    sub_mat = []

    # Convert substitution matrix into ints
    for row in sub_m:
        int_row = [int(numeric_string) for numeric_string in row]
        sub_mat.append(int_row)

    n = len(sequence1)  
    m = len(sequence2)


    # Generate matrix of zeros to store scores
    score_mat = np.zeros((m+1, n+1))

    # Generate the matrix for the paths
    path_array = np.zeros(score_mat.shape)
   
    # Calculate score table

    # first column
    for i in range(0, m + 1):
        path_array[i][0] = 2
        score_mat[i][0] = gap_penalty * i
    
    # first row
    for j in range(0, n + 1):
        path_array[0][j] = 1
        score_mat[0][j] = gap_penalty * j
    
    # Calulate the central values in the score mat
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if sequence1[j-1] == '-' or sequence2[i-1] == '-':
                match_score =  gap_penalty
            else:
                if sequence1[j-1] == 'A':
                    a = 0
                if sequence2[i-1] == 'A':
                    b = 0
                if sequence1[j-1] == 'C':
                    a = 1
                if sequence2[i-1] == 'C':
                    b = 1
                if sequence1[j-1] == 'T':
                    a = 2
                if sequence2[i-1] == 'T':
                    b = 2
                if sequence1[j-1] == 'G':
                    a = 3
                if sequence2[i-1] == 'G':
                    b = 3
                if a<4 and b<4:
                    match_score = sub_mat[a][b]
                else:
                    match_score = -100
            # # Calculate the score by checking the top, left, and diagonal cells
            # left_score = score_mat[i, j-1] + gap_penalty
            # up_score = score_mat[i-1, j] + gap_penalty
            # diag_score = score_mat[i-1, j-1] + gap_penalty
            match = score_mat[i - 1][j - 1] + match_score
            delete = score_mat[i - 1][j] + gap_penalty
            insert = score_mat[i][j - 1] + gap_penalty
            # Record the maximum score from the three possible scores calculated above
            cell_score = max(match, delete, insert)
            score_mat[i][j] = cell_score
            if cell_score == match:
                path_array[i][j] = 3
            if cell_score == delete:
                path_array[i][j] = 2
            if cell_score == insert:
                path_array[i][j] = 1
    
    
    # Initializes path array with 1 = left,  2 = up,  3 = diagonal

    path_array[0][0] = -1
    #print(path_array)
    # Use path array in order to find the optimal allignment
    row = m
    col = n
    trace = path_array[row][col]
    while trace != -1:
        trace = path_array[row][col]
        if(trace == 3):
            aligned_2 = sequence2[row-1] + aligned_2
            aligned_1 = sequence1[col-1] + aligned_1
            row -= 1
            col -= 1
        elif(trace ==2):
            aligned_2 = sequence2[row-1] + aligned_2
            aligned_1 = "-" + aligned_1
            row -= 1
        elif(trace ==1):
            aligned_2 = "-" + aligned_2
            aligned_1 = sequence1[col-1] + aligned_1
            col -= 1
        elif trace == -1:
            break
        else:
            print("ERROR: wrong symbol")
            break
    score = score_mat[m][n]
    return score, aligned_1, aligned_2


# for i in range(1, 11):
#     filename = "align_data_"+str(i)+".txt"
#     fasta1 = "prob4_data\data_"+str(i)+"\seq1.fasta"
#     fasta2 = "prob4_data\data_"+str(i)+"\seq2.fasta"
#     matfile = "prob4_data\data_"+str(i)+"\substitution_matrix.txt"
#     w = open(filename, "w")
#     score, dna1, dna2 = optimal_global_alignment(fasta1, fasta2, matfile)
#     w.write(f"The optimal alignment score between given sequences has score {score} \n")
#     w.write(f"{dna1}\n")
#     w.write(f"{dna2}\n")

# w = open("align_data_i.txt", "a")
# score, dna1, dna2 = optimal_global_alignment("prob4_data\data_10\seq1.fasta", "prob4_data\data_10\seq2.fasta", "prob4_data\data_10\substitution_matrix.txt")
# print(f"{score}\n{dna1}\n{dna2}")
# w.write(f"The optimal alignment score between given sequences has score {score} \n")
# w.write(f"{dna1}\n")
# w.write(f"{dna2}\n")
