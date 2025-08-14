# -*- coding: utf-8 -*-
"""
Python 3.6
1 December 2021;; jonathan.lee@pfizer.com
@author: LEEJ322

python align_to_refs.py consensus_fasta reference_seqs_fasta
"""
import os
import csv
import sys
from subprocess import call, DEVNULL, STDOUT

from Bio import Align, pairwise2, Seq, SeqIO
import pandas as pd



class VariantSimilarity:

    def __init__(self, query_file, path, passQC=True, ref_file='ospA_protein.fas'):
        self.passQC = passQC
        self.ref_file = ref_file  # .fasta format
        self.query_file = query_file  # .fasta format
        self.sample = self.query_file.split('_prot')[0]
        self.query = None
        self.refs = {}
        self.sero_tbl = {}
        self.variant = None
        self.serotype = None
        self.allele_scores = None
        self.path = path
        self.run()

    def load_query(self):
        # load query consensus sequence (protein sequence)
        inpath = '{0}/aln/{1}/{2}'.format(self.path, self.sample, self.query_file)
        print('align_to_refs:load_query - inpath: ' + inpath)

        with open(inpath, 'r') as inf:
            for record in SeqIO.parse(inf, 'fasta'):
                self.query = str(record.seq)
                if self.query[-1] == '*':  # removes stop codon
                    self.query = self.query[:-1]
                print('align_to_refs:load_query() - query: ', self.query)

    def load_refs(self):
        with open('variants/' + self.ref_file, 'r') as inf: #ref file is ospA_protein.fas
            for record in SeqIO.parse(inf, 'fasta'):
                seq = str(record.seq)
                id = int(str(record.description).split('_')[-1])
                self.refs[seq] = id

    def import_serotypes(self):
        with open('variants/serotypes.csv', 'r') as inf:
            reader = csv.reader(inf)
            next(reader)
            for row in reader:
                self.sero_tbl[int(row[0])] = (row[1], row[2])

    def calc_identity(self, query, seq):
        total = 0
        for i in range(len(query)):
            if query[i] == seq[i]: total += 1
        return total / len(query)

    def write_references(self):
        theFile = '{0}/aln/{1}/{1}_variants_ref_file'.format(self.path, self.sample)
        print('align_to_refs:write_references() - ' + theFile)

        with open(theFile, 'w') as outf:
            for seq, id in self.refs.items():
                outf.write('>ospA_protein_{}\n'.format(id))
                outf.write(seq + '\n')

    def run_mafft(self):
        self.write_references()
        theFile = '{0}/aln/{1}/{1}_variants_ref_file'.format(self.path, self.sample)

        build_cmd = 'cat {3}/aln/{0}/{1} {4} > {3}/aln/{0}/{0}_mafft_input.fasta'.format(self.sample, self.query_file,
                                                                                         self.ref_file, self.path,
                                                                                         theFile)
        call(build_cmd, stdout=DEVNULL, stderr=STDOUT, shell=True)

        print('align_to_refs:run_mafft() - ' + build_cmd)
        mafft_cmd = 'mafft {0}/aln/{1}/{1}_mafft_input.fasta > {0}/aln/{1}/{1}_mafft_out.fasta'.format(self.path,
                                                                                                       self.sample)

        print('align_to_refs:run_mafft() - ' + mafft_cmd)

        call(mafft_cmd, stdout=DEVNULL, stderr=STDOUT, shell=True)

    # alignments using Mafft
    def align_seqs(self):
        # run mafft
        self.run_mafft()
        # run comparisons
        query = None
        self.allele_scores = {}
        infile = '{0}/aln/{1}/{1}_mafft_out.fasta'.format(self.path, self.sample)
        for record in SeqIO.parse(infile, 'fasta'):
            seq = str(record.seq)
            if not query:
                query = seq
            else:
                ref = int(str(record.id).split('_')[-1])
                identity = self.calc_identity(query, seq)
                self.allele_scores[ref] = identity

        outfile = '{0}/aln/{1}/{1}_alignments.csv'.format(self.path, self.sample)
        with open(outfile, 'w') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['variant_id', 'identity'])
            for k, v in self.allele_scores.items():
                writer.writerow([k, v])

    def write_report(self):
        report_file = '{0}/aln/{1}_report.csv'.format(self.path, self.sample)
        with open(report_file, 'w') as outf:
            writer = csv.writer(outf)
            writer.writerow(['variant', self.variant])
            if self.allele_scores:
                writer.writerow(['similarity', max(self.allele_scores.values())])
            writer.writerow(['species', self.species])
            writer.writerow(['serotype', self.serotype])
        print('analysis complete. results written to ' + report_file)

    def run(self):
        if self.passQC:
            print('align_to_refs:run() - comparing consensus to known sequences...')
            self.load_query()
            self.load_refs()
            self.import_serotypes()
            # check if already in set
            self.variant = self.refs.get(self.query, 0)
            if self.variant != 0:
                print('align_to_refs:run() - sequence matches variant {}'.format(self.variant))
                self.serotype = self.sero_tbl[self.variant][0]
                self.species = self.sero_tbl[self.variant][1]
                
                # produce new fasta file of the OspA protein sequence with predicted epitope residues capitalized
                epitope_file = '{0}/aln/{1}/{1}_prot_epitope.fasta'.format(self.path, self.sample)
                seq_epitope_df = pd.read_csv('variants/ospA_protein_withEpitopes.csv')
                variant_id = 'OspA_protein_' + str(self.variant)
                mask = seq_epitope_df['id'] == variant_id
                # select the sequence corresponding to this variant
                capitalized_seq = seq_epitope_df.loc[mask, 'sequence'].iloc[0]
                with open(epitope_file, 'w') as epi:
                    epi.write('>consensus_with_epitope\n')
                    epi.write(capitalized_seq + '\n')

                # output csv file of the predicted epitope residues for this variant, along with RSA values of each residue
                epitope_csv = '{0}/aln/{1}/{1}_prot_epitope.csv'.format(self.path, self.sample)
                rsa_epitope_df = pd.read_csv('variants/bepipred_predicted_epitopes_rsa.csv')
                variant_df = rsa_epitope_df[rsa_epitope_df['protein'] == variant_id]
                variant_df.to_csv(epitope_csv,index=False)

                # output csv file of the vaccine variant and serotype that best matches the input variant
                epitope_analysis_csv = '{0}/aln/{1}/{1}_prot_epitope_analysis.csv'.format(self.path, self.sample)
                epi_analysis_df = pd.read_csv('variants/epitope_analysis.csv')
                epi_variant_df = epi_analysis_df[epi_analysis_df['protein'] == self.variant]
                row_dict = epi_variant_df.iloc[0].to_dict()
                with open(epitope_analysis_csv,'w',newline='') as f:
                    writer = csv.writer(f)
                    for key, value in row_dict.items():
                        writer.writerow([key,value])

                # output pymol session of the variant that the input genome corresponds to
                #pymol_file = "OspA_protein_" + str(self.variant) + ".pse"

            else:
                print('align_to_refs:run() - no match to known variants...')
                print('align_to_refs:run() - aligning to known sequences...')
                self.align_seqs()
                self.variant = max(self.allele_scores, key=self.allele_scores.get)
                self.serotype = self.sero_tbl[self.variant][0]
                self.species = self.sero_tbl[self.variant][1]
                print('align_to_refs:run() - highest similarity to {}'.format(self.variant))
        else:
            print('align_to_refs:run() - insufficient reads aligned or incomplete ospA sequence')
            self.variant = 'Unknown'
            self.serotype = 'Unknown'
            self.species = 'Unknown'
            # bad sequencing, not able to reconstruct full gene sequence during alignment and assembly
        self.write_report()


def main():
    VariantSimilarity(sys.argv[1])


if __name__ == '__main__':
    main()
