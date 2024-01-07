# importing the necessary modules
from Bio import SeqIO
import os
import re
import sys
import random

gene_target_file = open(sys.argv[1])  # read the contents of a list of genes from the command line
gene_target_list = gene_target_file.read().split('\n')  # split the contents of the file into a list based on a newline character
motif_file = open(sys.argv[2])  # read the contents of a list of motifs from the command line
motif_list = motif_file.read().split('\n')  # split the contents of the file into a list based on a newline character


def gene_list():
    """
    A function that iterates through each line of all gff3 files in the
    current working directory and filters for 'gene' sequence features.

    :return: list of all the rows in all the in .gff3 files in the current
    working directory that have "gene" as the sequence feature
    """
    gene = []
    for name in os.listdir("."):
        if name.endswith(".gff3"):
            annotation_file = open(name)
            for i in annotation_file.readlines():
                if i[0] != '#':
                    dataline_annotation = i.split('\t')
                    if dataline_annotation[2] == "gene":
                        gene.append(dataline_annotation)
    return (gene)


gene_list = gene_list()


def gene_target_hits():
    """
    A function that iterates through gene_list and gene_target_list and
    keeps all elements in gene_list where the gene id matches one of
    the gene id's in gene_target_list

    :return: list of genes from gene_list that matches with
    gene_target_list
    """
    gene_target_hits = []
    for i in gene_list:
        id = i[8]  # element that contains the id for the iterating gene
        match = re.search(r'gene:([^;]+)',id)  # regular expression to match the specific gene id for the iterating element in gene_list
        if match:
            gene_id = match.group(1)
            for j in gene_target_list:
                if gene_id == j:
                    gene_target_hits.append(i)  # append the iterating element in gene_list to gene_target_list if gene_id matches a iterating element in gene_target_list
    return (gene_target_hits)


gene_list_filtered_for_matches_with_gene_target_list = gene_target_hits()


def sequence_file():
    """
    A function that creates a new file called outfile_sequences.fa and
    writes the contents of all fasta files in the current working
    directory into outfile_sequences.fa


    :return:
    """
    out = open("outfile_sequences.fa", "w")
    for name in os.listdir("."):
        if name.endswith(".fa"):
            sequence_file = open(name)
            out.write(sequence_file.read().rstrip() + '\n')
            sequence_file.close()
    out.close()


sequence_file()



def chromosome_number_sequence():
    """
    A function that creates a dictionary called 'chromosomeID_seq' where the key is the
    chromosome number from the fasta sequences in outfile_sequences.fa and the value is the
    nucleotide sequence belonging to the specific chromosome

    :return: A dictionary where the key corresponds to a chromosome number and the value
    is the chromosome sequence
    """
    sequences = SeqIO.parse('outfile_sequences.fa', 'fasta')
    chromosomeID_seq = {}
    for i in sequences:
        chromosomeID_seq[
            int(i.id)] = i.seq  # i.id refers to the chromosome number and i.seq is the nucelotide sequence of that chromosome
    return (chromosomeID_seq)


chromosome_number_sequence = chromosome_number_sequence()


def random_genes():
    """
    A function that randomly chooses 593 genes from gene_list that are not in
    gene_list_filtered_for_matches_with_gene_target_list.
    :return: A list of 593 randomly chosen genes from gene_list
    """
    random_gene_list = []
    for i in gene_list_filtered_for_matches_with_gene_target_list:
        for j in gene_list:
            if i != j:
                random_gene_list.append(j)
    random_gene = random.sample(random_gene_list, 593)  # using the random module to create a list of 593 random genes from random_gene_list
    return (random_gene)


random_genes1 = random_genes()
random_genes2 = random_genes()
random_genes3 = random_genes()
random_genes4 = random_genes()
random_genes5 = random_genes()


def promoter_region(gff3_gene_file):
    """
    A function that creates a dictionary where the key is a gene ID
    for the given iterating element from a list of genes from a gff3 file
    and the value is the promoter region for that specific gene.
    :param gff3_gene_file:
    :return: A dictionary where the key is a gene ID and the value is the upstream
    promoter region for the gene
    """
    gene_promoter_region = {}
    for i in gff3_gene_file:
        chromosome_for_gene = int(i[0])  # chromosome number for the iterating element in gff3_gene_file
        start = int(i[3])  # start site for the iterating element which represents a gene
        end = int(i[4])  # end site for the iterating element which represents a gene
        sequence_for_gene = chromosome_number_sequence[chromosome_for_gene]  # using the chromosome_number_sequence dictionary to assign the specific chromosome sequence based on the iterating element to sequence_for_gene
        id = i[8]
        match = re.search(r'gene:([^;]+)', id)  # getting the gene id for the iterating element
        if match:
            gene_id = match.group(1)  # assign the gene id to gene_id
            if i[6] == "+":  # if the iterating gene is on the positive strand
                gene_sequence = sequence_for_gene[start - 500: start]  # index the promoter region upstream of the start and assign to gene_sequence
                reverse_gene_sequence = gene_sequence[::-1]  # reverse gene_sequence
                gaps = re.search(r'N{100}',str(reverse_gene_sequence))  # match the first instance of 100 consecutive N's which represents a gap
                if gaps:  # if a gap is found in the iterative promoter region
                    gene_sequence_gaps = gene_sequence[len(gene_sequence) - gaps.start():]  # index gene_sequnce to first upstream gap
                    gene_promoter_region[gene_id] = gene_sequence_gaps  # create a new key-value pair with the iterating gene id as the key and gene_sequence_gaps as the value
                else:  # if a gap is not found in the iterating promoter region
                    gene_promoter_region[gene_id] = gene_sequence  # create a new key-value pair for gene_promoter_region with the iterating gene id as the key and the gene_sequence as the value
            elif i[6] == "-":  # if the iterating gene is on the negative strand
                gene_sequence = sequence_for_gene[end:end + 500]  # indexing the upstream promoter region for the negative strand
                rc_seq = gene_sequence.reverse_complement()  # take the reverse complement of gene_sequence and assign to rc_seq
                reverse_gene_sequence = rc_seq[::-1]  # reverse rc_seq
                gaps = re.search(r'N{100}',str(reverse_gene_sequence))  # match the first instance of 100 consecutive N's which represents a gap
                if gaps:  # if a gap is found in the iterative promoter region
                    gene_sequence_gaps = rc_seq[len(rc_seq) - gaps.start():]  # index gene sequence to the first upstream gap of rc_seq
                    gene_promoter_region[gene_id] = gene_sequence_gaps  # create a new key value pair for gene_promoter_region where the iterating gene_id is the key and the gene_sequence_gaps as the value
                else:  # if a gap is not found in the iterative promoter region
                    gene_promoter_region[gene_id] = rc_seq  # create a new key-value pair for gene_promoter_region with the iterating gene id as the key and the rc_seq as the value
    return (gene_promoter_region)


selected_run = promoter_region(gene_list_filtered_for_matches_with_gene_target_list)

random1 = promoter_region(random_genes1)
random2 = promoter_region(random_genes2)
random3 = promoter_region(random_genes3)
random4 = promoter_region(random_genes4)
random5 = promoter_region(random_genes5)


def promoter_motif_match(geneID_with_promoter_region):
    """
    A function that formats the motifs from motif_list and
    creates a dictionary where the key is the motif and the value is
    the number of times the motif is matched in each of the
    promoter regions for a given input of genes
    :param geneID_with_promoter_region:
    :return: A dictionary where the key is a motif from motif_list
    and the value is the count of the occurrences of the motif in the
    promoter regions from geneID_with_promoter_region
    """
    formatted_motifs = [] #list of motifs from motif_list that are capitalized
    for motifs in motif_list:
        formatted_motifs.append(motifs.upper())

    motif_matches = {}
    for motifs in formatted_motifs: #initiaing a dictionary where each key is a motif from formatted_motifs and the value is zero
        motif_matches[motifs] = 0
    for promoter_region in geneID_with_promoter_region.values(): # count each instance of every motif in formatted_motifs for each promoter sequence in geneID_with_promoter_region
        for matching_motifs in formatted_motifs:
            search = re.findall(matching_motifs, str(promoter_region))
            if search:
                for i in search: #for each match increment the value for the current matching_motifs by 1
                    motif_matches[matching_motifs] += 1

    return (motif_matches)


final_selected = promoter_motif_match(selected_run)
final_random1 = promoter_motif_match(random1)
final_random2 = promoter_motif_match(random2)
final_random3 = promoter_motif_match(random3)
final_random4 = promoter_motif_match(random4)
final_random5 = promoter_motif_match(random5)


def output_file():
    """
    A function that creates a combined dictionary from
    final_selected, final_random1, final_random2, final_random3,
    final_random4, and final_random5. The keys are the motifs
    and the first value is the counts of the motif found in the
    upstream promoter region of the selected genes and
    the next 5 values are the counts for 5 instances of 593
    randomly selected genes
    :return:
    """
    combined_dict = {}

    for key in final_selected: # looping through each key in final_selected and creating a new dictionary where each key is a motif and the values are the count of the occurrences of that motif in the upstream promoter region of the selected genes and the 5 randomly chosen genes. Each value is separated by a tab
        combined_dict[key] = (str(final_selected[key]) + '\t' + str(final_random1[key]) + '\t' + str(final_random2[key]) + '\t' + str(final_random3[key]) + '\t' + str(final_random4[key]) + '\t' + str(final_random5[key]))

    final_output_file_path = "final_output.txt"
    with open(final_output_file_path, 'w') as file: #creatng a file called final_output.txt and writing each motif and the counts to the file
        for motif, count in combined_dict.items():
            file.write(f"{motif}\t{count}\n")


output_file()
