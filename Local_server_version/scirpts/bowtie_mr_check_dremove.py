#!/usr/bin/env python
# bowtie_mr_check_dremove.py
# 2016-10-20

"""
bowtie_mr_check_dremove.py   [ --single | --paired ]  [ --DNase ]
                             --exp-ss=100  --final-ss=100
                             --read-len=36  [ --not-remove-duplicates ]
                             --stats=NA.dup.stats
                             Input.mr  > Output.dremove.bed

[ --single | --paired ]      Specify whether the reads are single-
                             or paired-end.
                             Default: single-end

--DNase                      Specify whether the reads are from a
                             DNase-seq experiment.
                             Default: OFF

--exp-ss=<int>               1/2 of estimated DNA fragment size in
                             the experiment for library preparation.
                             Used to determine the exact binding sites.
                             Ignored when --paired and/or --DNase is set.
                             Must be a non-negative integer.
                             Default: 100

--final-ss=<int>             Shift size to the 3' direction to reach the
                             exact binding sites in the OUTPUT BED file.
                             Must be a non-negative integer.
                             Default: 100

--read-len=<int>             Uniform read length in the OUTPUT BED file.
                             Must be a positive integer.
                             Default: 36

--not-remove-duplicates      Not to remove any duplicate reads or read pairs.
                             Default: OFF

--stats=<str>                File to write summary statistics about
                             removing duplicates to.
                             Default: NA.dup.stats

<Input.mr>                   Map results from Bowtie 1 (using default
                             format). The first 5 columns are required.
                             Use - to specify the standard input stream.
                             Default: -

<Output.dremove.bed>         Duplicates removed read positions will be
                             written to the standard output stream.

Check items:
1. Only one alignment for each read or read pair is allowed.
2. For a read pair, '+' record occurs before '-' record.
3. Negative coordinates won't be checked and might be written to the stdout.

Figure out the exact binding sites:
1. If --DNase is set, 5' end of each read (mate) is the exact binding site. Thus,
   two lines will be written for each read pair.
2. Else if --paired is set, the center point of deduced DNA fragment is the exact
   binding site. Thus, only one line is written for each read pair.
3. Else, the exact binding site of each read is its 5' end shifted <--exp-ss> bases
   downstream.

Strand information in <Output.dremove.bed>:
1. If --paired is set and --DNase is not set, the strand field in the output for each
   read pair is the same as that of mate 1.
2. Otherwise, the strand field in the output is the same as that of the read in question.

"""


from sys import argv, stderr, stdin, stdout
from getopt import getopt


def single_end(f, exp_ss=None, DNase=False, dup_remove=True):
    if not DNase and exp_ss is None:
        raise Exception("Experimental shift size must be provided when the case is not DNase-seq!")
    
    # A read position is represented by (chr, strand, 5' end position).
    pos = set()
    names = set()
    for line in f:
        items = line.rstrip("\r\n").split("\t")
        if items[0] in names:
            raise Exception("Multiple alignment records exist for the read %r" % items[0])
        names.add(items[0])
        if items[1] != "+" and items[1] != "-":
            raise Exception("Invalid strand field in the alignment record for %r: %s" % (items[0], items[1]))
        
        end_5 = int(items[3])
        if items[1] == "-":
            end_5 += len(items[4]) - 1
        site = (items[2], items[1], end_5)
        if site in pos and dup_remove:
            continue
        pos.add(site)
        
        if DNase:
            yield site
        else:
            if items[1] == "+":
                yield (items[2], items[1], end_5 + exp_ss)
            else:
                yield (items[2], items[1], end_5 - exp_ss)
    
    yield (len(pos), len(names))


def paired_end(f, DNase=False, dup_remove=True):
    # A read pair position is represented by (chr, 5 end position in + strand, 5 end position in - strand)
    pos = set()
    names = set()
    unpaired = {}
    for line in f:
        items = line.rstrip("\r\n").split("\t")
        if not (items[0].endswith("/1") or items[0].endswith("/2")):
            raise Exception("Invalid name for a mate of a read pair: %r" % items[0])
        name = items[0][:-2]
        if name in names:
            raise Exception("Multiple alignment records exist for the read pair %r" % name)
        if name not in unpaired:
            if items[1] != "+":
                raise Exception("The 1st alignment record for the read pair %r is not in + strand!" % name)
            unpaired[name] = (items[0][-1], items[2], int(items[3]))
        else:
            if items[1] != "-":
                raise Exception("The 2nd alignment record for the read pair %r is not in - strand!" % name)
            mate, ch, plus_5_end = unpaired.pop(name)
            if mate == items[0][-1]:
                raise Exception("The names of the two mates of read pair %r are not consistent to each other!" % name)
            if ch != items[2]:
                raise Exception("The chromatin fields in the two alignment records for read pair %r are not identical!" % name)
            minus_5_end = int(items[3]) + len(items[4]) - 1
            if plus_5_end > minus_5_end:
                raise Exception("The orientation of the two mates of read pair %r are not correct!" % name)
            
            names.add(name)
            site = (ch, plus_5_end, minus_5_end)
            if site in pos and dup_remove:
                continue
            pos.add(site)
            
            if DNase:
                yield (ch, "+", plus_5_end)
                yield (ch, "-", minus_5_end)
            else:
                strand = "+" if mate == "1" else "-"
                yield (ch, strand, (plus_5_end + minus_5_end) // 2)
    
    if len(unpaired) != 0:
        i, j = unpaired.items()[0]
        raise Exception("Some alignment records remain unpaired!\nInstance: the record for the mate %r" % (i + "/" + j[0]))
    
    yield (len(pos), len(names))


def util_output(hit, final_ss, read_len, cnt):
    ch, strand, pos = hit
    if strand == "+":
        start = pos - final_ss
        end = start + read_len
    elif strand == "-":
        end = pos + final_ss + 1
        start = end - read_len
    stdout.write("%s\t%d\t%d\tHIT%d\t0\t%s\n" % (ch, start, end, cnt, strand))


def main():
    mode = "single"
    DNase = False
    exp_ss = 100
    final_ss = 100
    read_len = 36
    stats = "NA.dup.stats"
    mr = "-"
    dup_remove = True
    try:
        opts, args = getopt(argv[1:], "h", ["single", "paired", "DNase", "exp-ss=", "final-ss=",
                                            "read-len=", "not-remove-duplicates", "stats=", "help"])
        if len(args) != 0:
            if len(args) > 1:
                raise Exception("Multiple input files specified: %r" % args)
            mr = args[0]
        for i, j in opts:
            if i == "-h" or i == "--help":
                stdout.write(__doc__)
                exit(0)
            elif i == "--single":
                mode = "single"
            elif i == "--paired":
                mode = "paired"
            elif i == "--DNase":
                DNase = True
            elif i == "--exp-ss":
                try:
                    exp_ss = int(j)
                    assert exp_ss >= 0
                except (ValueError, AssertionError):
                    raise Exception("Invalid argument for %s: %s" % (i, j))
            elif i == "--final-ss":
                try:
                    final_ss = int(j)
                    assert final_ss >= 0
                except (ValueError, AssertionError):
                    raise Exception("Invalid argument for %s: %s" % (i, j))
            elif i == "--read-len":
                try:
                    read_len = int(j)
                    assert read_len > 0
                except (ValueError, AssertionError):
                    raise Exception("Invalid argument for %s: %s" % (i, j))
            elif i == "--not-remove-duplicates":
                dup_remove = False
            elif i == "--stats":
                stats = j
            else:
                raise Exception("Internal errors occur when parsing command line arguments.")
    except Exception as e:
        stderr.write("%s\n" % e)
        stderr.write("Try 'bowtie_mr_check_dremove.py --help' for more information.\n")
        exit(1)
    
    f = stdin if mr == "-" else open(mr)
    with f:
        cnt = 1
        if mode == "single":
            hits = single_end(f, exp_ss, DNase, dup_remove)
        elif mode == "paired":
            hits = paired_end(f, DNase, dup_remove)
        for hit in hits:
            if len(hit) == 3:
                util_output(hit, final_ss, read_len, cnt)
                cnt += 1
            else:
                uniq, total = hit
    
    dup = total - uniq
    unit = "reads" if mode == "single" else "read pairs"
    with open(stats, "w") as w:
        w.write("Total %s: %d\n" % (unit, total))
        w.write("Unique %s: %d (%.2f%s)\n" % (unit, uniq, uniq * 100.0 / total, "%"))
        w.write("Duplicate %s: %d (%.2f%s)\n" % (unit, dup, dup * 100.0 / total, "%"))





if __name__ == '__main__':
    main()

