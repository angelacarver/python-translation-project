#! /usr/bin/env python3

import sys

def pop_next_codon(sequence):
    """Removes and returns the first 3 bases.

    Returns a tuple of a string of the first three bases and a string of the remaing sequence.
    """
    codon = sequence[0:3]  #takes the first three bases
    remaining_seq = sequence[3:]   #the rest of the sequence
    return codon, remaining_seq   #returns the two parts of the sequence

joinedgenes = ""
def translate_sequence(rna_sequence, genetic_code):
    """Translates a sequence of RNA into a sequence of amino acids.

    Translates `rna_sequence` into string of amino acids, according to the
    `genetic_code` given as a dict. Translation begins at the first position of
    the `rna_sequence` and continues until the first stop codon is encountered
    or the end of `rna_sequence` is reached.

    If `rna_sequence` is less than 3 bases long, or starts with a stop codon,
    an empty string is returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    str
        A string of the translated amino acids.
    """
    rna_sequence = rna_sequence.upper() #makes it all upper case
    amino_acid_list = []  #setting up variable
    while True:
        if len(rna_sequence) <3:
            break #if the sequence is less than three bases long, returns empty string
        codon, remaining_seq = pop_next_codon(rna_sequence) #this function pop_next_codon is defined earlier
        rna_sequence = remaining_seq  #this is the sequence after the first three bases as seen in the function pop_next_codon
        aa = genetic_code[codon]
        if aa == "*":
            break  #if the sequence contains stop codon, returns empty string
        amino_acid_list.append(aa)
    return "".join(amino_acid_list)

def get_all_translations(rna_sequence, genetic_code):
    """Get a list of all amino acid sequences encoded by an RNA sequence.

    All three reading frames of `rna_sequence` are scanned from 'left' to
    'right', and the generation of a sequence of amino acids is started
    whenever the start codon 'AUG' is found. The `rna_sequence` is assumed to
    be in the correct orientation (i.e., no reverse and/or complement of the
    sequence is explored).

    The function returns a list of all possible amino acid sequences that
    are encoded by `rna_sequence`.

    If no amino acids can be translated from `rna_sequence`, an empty list is
    returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    list
        A list of strings; each string is an sequence of amino acids encoded by
        `rna_sequence`.
    """
    rna_sequence = rna_sequence.upper() #makes the sequence all capital letters
    number_of_bases = len(rna_sequence) #gives the length of the rna sequence
    last_codon_index = number_of_bases - 3 #shows where the last codon in rna_sequence starts
    if last_codon_index < 0:   #tests whether rna_sequence is long enough to contain any codons
        return[]  #if it's too short, return empty list
    amino_acid_seq_list = []
    for base_index in range(last_codon_index +1):  #gives number of bases to loop through, not sure why last_codon_index + 1
        codon = rna_sequence[base_index: base_index +3] #indicates that each codon is three bases starting at the base_index
        if codon == "AUG": #if start codon appears
            aa_seq = translate_sequence(         #use the translation function we defined to translate sequence
                rna_sequence = rna_sequence[base_index:],
                genetic_code = genetic_code)
            if aa_seq: #if aa_seq happens
                amino_acid_seq_list.append(aa_seq) #add this to the list we defined earlier
    return amino_acid_seq_list  #return the list of codons


reverse = "string"
reverse_upper = "string"
def get_reverse(sequence):
    """Reverse orientation of `sequence`.

    Returns a string with `sequence` in the reverse order.

    If `sequence` is empty, an empty string is returned.

    Examples
    --------
    >>> get_reverse('AUGC')
    'CGUA'
    """
    reverse = sequence[::-1]
    reverse_upper = (reverse.upper())
    return reverse_upper

def get_complement(sequence):
    """Get the complement of a `sequence` of nucleotides.

    Returns a string with the complementary sequence of `sequence`.

    If `sequence` is empty, an empty string is returned.

    Examples
    --------
    >>> get_complement('AUGC')
    'UACG'
    """
    complementary_strand = ""
    for base in sequence :
        if base == "A" :
            complementary_strand += "U"
        
        elif base == "a" :
            complementary_strand += "U"

        elif base == "U" :
            complementary_strand += "A"

        elif base == "u" :
            complementary_strand += "A"

        elif base == "G" :
            complementary_strand += "C"

        elif base == "g" :
            complementary_strand +="C"

        elif base == "C" :
            complementary_strand += "G"

        elif base == "c" :
            complementary_strand += "G"

        else :
            print("Wrong input")
            break

    return complementary_strand

reverse_2 = "string"
compstrand = ""
def reverse_and_complement(sequence):
    """Get the reversed and complemented form of a `sequence` of nucleotides.

    Returns a string that is the reversed and complemented sequence
    of `sequence`.

    If `sequence` is empty, an empty string is returned.

    Examples
    --------
    >>> reverse_and_complement('AUGC')
    'GCAU'
    """
    reve = get_reverse(sequence)
    reve_comp = get_complement(reve)
    return reve_comp


def get_longest_peptide(rna_sequence, genetic_code):
    """Get the longest peptide encoded by an RNA sequence.

    Explore six reading frames of `rna_sequence` (the three reading frames of
    `rna_sequence`, and the three reading frames of the reverse and complement
    of `rna_sequence`) and return (as a string) the longest sequence of amino
    acids that it encodes, according to the `genetic_code`.

    If no amino acids can be translated from `rna_sequence` nor its reverse and
    complement, an empty string is returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    str
        A string of the longest sequence of amino acids encoded by
        `rna_sequence`.
    """
    peptides = get_all_translations(rna_sequence = rna_sequence,
        genetic_code = genetic_code) #use the get_all_translations function to translate the base rna_sequence
    rev_comp_seq = reverse_and_complement(rna_sequence) #use the reverse_and_complement function that we defined to get the reverse/complement of rna_sequence
    rev_comp_peptides = get_all_translations(rna_sequence = rev_comp_seq,
        genetic_code = genetic_code) #use the get_all_translations function to translate the reverse/complement
    peptides += rev_comp_peptides #add this translation to peptides variable
    if not peptides:
        return ""    #unsure what this section means
    if len(peptides) <2:  #if doesn't contain a two codons, return the first amino acid
        return peptides[0]
    most_number_of_bases = -1
    longest_peptide_index = -1
    for peptide_index, aa_seq in enumerate(peptides):  #enumerate function tracks iterations of the loop
        if len(aa_seq) > most_number_of_bases: #if more than -1 bases?
            longest_peptide_index = peptide_index  #the counter is equal to the longest peptide index
            most_number_of_bases = len(aa_seq)  #change the most number of bases so that it moves to the next index next time
    return peptides[longest_peptide_index] #return the last peptide

if __name__ == '__main__':
    genetic_code = {'GUC': 'V', 'ACC': 'T', 'GUA': 'V', 'GUG': 'V', 'ACU': 'T', 'AAC': 'N', 'CCU': 'P', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'CAU': 'H', 'AAU': 'N', 'AGU': 'S', 'GUU': 'V', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'CCA': 'P', 'ACA': 'T', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'UGC': 'C', 'CAG': 'Q', 'GAU': 'D', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': '*', 'GGA': 'G', 'UAC': 'Y', 'GAC': 'D', 'UAG': '*', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'GGC': 'G', 'AUG': 'M', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'GCU': 'A', 'GAA': 'E', 'AUU': 'I', 'UUG': 'L', 'UUA': 'L', 'UGA': '*', 'UUC': 'F'}
    rna_seq = ("AUG"
            "UAC"
            "UGG"
            "CAC"
            "GCU"
            "ACU"
            "GCU"
            "CCA"
            "UAU"
            "ACU"
            "CAC"
            "CAG"
            "AAU"
            "AUC"
            "AGU"
            "ACA"
            "GCG")
    longest_peptide = get_longest_peptide(rna_sequence = rna_seq,
            genetic_code = genetic_code)
    assert isinstance(longest_peptide, str), "Oops: the longest peptide is {0}, not a string".format(longest_peptide)
    message = "The longest peptide encoded by\n\t'{0}'\nis\n\t'{1}'\n".format(
            rna_seq,
            longest_peptide)
    sys.stdout.write(message)
    if longest_peptide == "MYWHATAPYTHQNISTA":
        sys.stdout.write("Indeed.\n")
