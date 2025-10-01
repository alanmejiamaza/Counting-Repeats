#!/usr/bin/env python

import argparse as ap
from typing import IO


def reverse_complement(seq, rna=False):
    complements = {'A': 'T',
                   'T': 'A',
                   'G': 'C',
                   'C': 'G',
                   'N': 'N'}

    if rna:
        seq.replace("U", "T")
    # Reverse complement
    rc = "".join(complements[s] for s in seq[::-1])
    if rna:
        rc.replace("T", "U")

    return rc


class Read:
    class StartMotifError(Exception):
        """Exception raised when the start motif of a subsequence cannot be found."""
        pass

    def __init__(self, name: str, seq: str, qual: str = None, rna=False):
        self.subseq = None
        self.name = str(name)
        self.seq = seq
        self.qual = qual
        self.rna = rna

    def __repr__(self):
        return f"{self.__class__}; {self.name}; {self.seq}; current subseq: {self.subseq}"

    def __str__(self):
        if self.subseq:
            return self.subseq
        else:
            return self.seq

    def __eq__(self, other):
        # if isinstance(other, Read):
        #     return str(self) == str(other)
        # else:
        #
        #     return False
        try:
            return str(self) == str(other)
        except ValueError:
            return False

    def __hash__(self):
        return hash(str(self))

    def subset(self, start, stop):
        # Find first instance of start motif. Raise MotifStartError if the start motif cannot be found.
        subset_start = self.seq.find(start)
        if subset_start == -1:
            raise self.StartMotifError(f"Could not find {start} in read {self.name}")
        # Find the stop motif downstream from the start motif. If not found, read until the end of the sequence
        subset_stop = self.seq.find(stop, subset_start)
        if subset_stop == -1:
            self.subseq = self.seq[subset_start:]
        else:
            self.subseq = self.seq[subset_start: subset_stop + len(stop)]

    def reorient(self, motif, subseq=False):
        if subseq and self.subseq:
            s = self.subseq
        else:
            s = self.seq

        motif_count = s.count(motif)
        motif_rc_count = s.count(reverse_complement(motif))

        if motif_rc_count > motif_count:
            self.seq = reverse_complement(self.seq, rna=self.rna)
            if self.subseq:
                self.subseq = reverse_complement(self.subseq, rna=self.rna)

    def fastx(self):
        fastx = [self.name, str(self), '+', self.qual]
        if self.qual:
            return fastx
        else:
            return fastx[:2]

    @classmethod
    def from_record(cls, rec: Union[list, tuple]):
        # Create instance from a flattened fasta/q record
        name, seq = rec[:2]
        if len(rec) == 4:
            qual = rec[3]
        else:
            qual = None

        return cls(name=name, seq=seq, qual=qual)


class Pair:
    def __init__(self, r1: Read, r2: Read):
        self.r1 = r1
        self.r2 = r2

    def __iter__(self):
        yield self.r1
        yield self.r2

    def __str__(self):
        return self.r1.name

    def is_consensus(self):
        return self.r1 == self.r2


class Aggregates(dict):
    table_header = ['sample', 'sequence', 'count', 'proc_freq', 'total_freq']

    def __init__(self, sample):
        self.sample = sample
        self.table = None
        super().__init__()

    def calculate_freq(self, count_unproc_reads: int = 0):
        if not self.table:
            self.to_table()
        # Calculate frequency amongst processed reads and total reads
        # self.table_header.append("proc_freq")
        # self.table_header.append("total_freq")
        count_proc_reads = sum(row[2] for row in self.table)
        for row in self.table:
            freq_p = row[2] / count_proc_reads
            freq_t = row[2] / (count_proc_reads + count_unproc_reads)
            row.append(freq_p)
            row.append(freq_t)

    def to_table(self):
        agg_table = [[self.sample, s, c['count']] for s, c in self.items()]
        agg_table = sorted(agg_table, key=lambda x: x[2], reverse=True)

        self.table = agg_table

    def to_file(self, handle: IO):
        handle.write("\t".join(self.table_header) + '\n')
        for row in self.table:
            handle.write("\t".join(str(x) for x in row) + '\n')

    def assign_reads(self):
        header = ['sample', 'seq', 'read_id']
        table = [header]
        for seq, stats in self.items():
            for r in stats['reads']:
                table.append([self.sample, seq, r])

        return table


class Unprocessed(list):
    table_header = ['sample', 'read_id', 'failure']

    def __init__(self, sample: str):
        self.sample = sample
        super().__init__()
        self.append(self.table_header)

    def __iter__(self):
        for row in super().__iter__():
            try:
                yield [self.sample] + row
            except TypeError:
                yield [self.sample] + [row]


def write_to_file(table: list, handle: IO):
    for row in table:
        handle.write("\t".join(str(r) for r in row) + '\n')


def preprocess_read(r: Read, start, stop, test_orient=None):
    if test_orient:
        r.reorient(test_orient, subseq=True)
    r.subset(start, stop)


def characterize(read: Read, agg: dict):
    try:
        agg[read.subseq]['count'] += 1
        agg[read.subseq]['reads'].append(read.name)
    except KeyError:
        agg[read.subseq] = {"count": 1, "reads": [read.name]}


def main(sample: str, fq1, start: str, stop: str, fq2=None, test_orient: str = None, parse_from_file=False):
    agg = Aggregates(sample)
    unproc = Unprocessed(sample)

    if fq2:
        for mate1, mate2 in zip(fq1, fq2):
            if parse_from_file:
                mate1 = mate1.strip().split('\t')
                mate2 = mate2.strip().split('\t')
            m1 = Read.from_record(mate1)
            m2 = Read.from_record(mate2)

            pair = Pair(m1, m2)

            for read in pair:
                try:
                    preprocess_read(read, start, stop, test_orient=test_orient)
                    # Only analyze the subsequence if both mates of the pair exactly match
                    if pair.is_consensus():
                        characterize(pair.r1, agg)
                    else:
                        unproc.append([pair, 'mismatch'])  # Cannot characterize - Mates' subseq differ
                except Read.StartMotifError:
                    unproc.append([pair, 'start_motif_not_found'])  # Cannot characterize - Motif not found

    else:
        for record in fq1:
            if parse_from_file:
                record = record.strip().split('\t')
            read = Read.from_record(record)
            try:
                preprocess_read(read, start, stop, test_orient=test_orient)
                characterize(read, agg)
            except Read.StartMotifError:
                unproc.append([read, 'start_motif_not_found'])

    agg.calculate_freq(len(unproc))

    return agg, unproc


if __name__ == "__main__":
    parser = ap.ArgumentParser(description="Aggregate (sub)sequences of reads")
    parser.add_argument("sample", help="Sample name being analyzed.")
    parser.add_argument("fq1", help="Flattened fasta/q file")
    parser.add_argument("start", help="Motif that is the start of the subsequence to characterize.")
    parser.add_argument("stop", help="Stop motif that terminates the desired subsequence.")
    parser.add_argument("--fq2", help="Balanced mates of first files reads.")
    parser.add_argument("--test_orient", help="Motif to orient read sequence. Checks for reverse complementing.")
    parser.add_argument("--output_prefix", default="", help="Additional strings for output file prefixes.")

    args = parser.parse_args()
    # Set up variables for single or paired sequencing
    seq1 = open(args.fq1, 'r')
    if args.fq2:
        seq2 = open(args.fq2, 'r')
        output_prefix = f"{args.sample}{args.output_prefix}_paired_consensus"
    else:
        seq2 = None
        output_prefix = f"{args.sample}{args.output_prefix}_single"

    aggregate, unprocessed = main(args.sample, seq1, args.start, args.stop, fq2=seq2, test_orient=args.test_orient,
                                  parse_from_file=True)
    seq1.close()
    if args.fq2:
        seq2.close()

    output_files = [f"{output_prefix}.aggregate_seq.txt",
                    f"{output_prefix}.unprocessed.txt",
                    f"{output_prefix}.processed_reads.txt"]
    output_data = [aggregate.table,
                   unprocessed,
                   aggregate.assign_reads()]

    for data, out_file in zip(output_data, output_files):
        with open(out_file, 'w') as fw:
            write_to_file(data, fw)
