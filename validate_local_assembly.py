import pybedtools
import pysam
import tempfile
import shutil
import subprocess
import os.path
import time
import sqlite3
import sys
import argparse

def read_sam_file(sam_file_name):
    samfile = pysam.Samfile(sam_file_name, "r")
    records_by_name = {}
    for record in samfile:
        record_name = record.qname
        if not record_name in records_by_name:
            records_by_name[record_name] = []
        records_by_name[record_name].append(record)
    return (samfile, records_by_name)

def get_dels_from_cigar(rname, cigar, pos):
    BAM_CSOFT_CLIP = 4
    BAM_CHARD_CLIP = 5
    BAM_CDEL = 2
    i = 0
    deletions_in_cigar = []
    while cigar[i][0] in [BAM_CSOFT_CLIP, BAM_CHARD_CLIP]:
        i = i + 1
    while i < len(cigar):
        if cigar[i][0] == BAM_CDEL:
            deletion_interval = pybedtools.Interval(rname, pos, pos + cigar[i][1])
            deletions_in_cigar.append(deletion_interval)
            pos = pos + cigar[i][1]
        else:
            pos = pos + cigar[i][1]
        i = i + 1
    return deletions_in_cigar
            

def left_hard_clipping(record):
    BAM_CHARD_CLIP = 5
    return 0 if record.cigar[0][0] != BAM_CHARD_CLIP else record.cigar[0][1]

def get_deletion_intervals(samfile, sam_records_by_read):
    BAM_CSOFT_CLIP = 4
    BAM_CHARD_CLIP = 5
    deletion_intervals = []
    for qname in sam_records_by_read:
        records = sorted(sam_records_by_read[qname], key=lambda record: (record.qend + left_hard_clipping(record)))
        records = sorted(records, key=lambda record: (record.qstart + left_hard_clipping(record)))
        records = sorted(records, key=lambda record: record.tid)
        print map(str, records)
        qname_idx = 1
        filtered_records = []

        for i in range(len(records)):
            print "record {0} qstart = {1} quend {2}".format(i, records[i].qstart + left_hard_clipping(records[i]), records[i].qend + left_hard_clipping(records[i]))
            if i == 0 or records[i].qstart + left_hard_clipping(records[i]) >= records[i-1].qend + left_hard_clipping(records[i-1]) - 50:
                print "keeping records {0}\n".format(records[i])
                filtered_records.append(records[i])
        
        records = filtered_records

        for i in range(len(records)):
            deletion_intervals = deletion_intervals + get_dels_from_cigar(samfile.getrname(records[i].tid), records[i].cigar, records[i].pos)
            is_reverse = True
            if i < (len(records) - 1):
                a = records[i]                
                is_reverse = a.is_reverse                
                b = records[i+1]  
                
                print qname
                print str(a.pos) + "-" + str(a.aend)
                print str(a.qstart) + "-" + str(a.qend)
                print a.cigar
                print str(b.pos) + "-" + str(b.aend)
                print str(b.qstart) + "-" + str(b.qend)
                print b.cigar          
                left_read_right_clipping = 0
                i = len(a.cigar)
                while i > 0 and a.cigar[i-1][0] in [BAM_CSOFT_CLIP, BAM_CHARD_CLIP]:
                    left_read_right_clipping += a.cigar[i-1][1]
                    i = i - 1
                print "left read clipping: {0}".format(left_read_right_clipping)
                right_read_left_clipping = 0
                i = 0
                while i < len(b.cigar) and b.cigar[i][0] in [BAM_CSOFT_CLIP, BAM_CHARD_CLIP]:
                    right_read_left_clipping += b.cigar[i][1]
                    i = i + 1
                print "right read clipping: {0}".format(right_read_left_clipping)
                if a.tid != b.tid:
                    continue
                if is_reverse != b.is_reverse:
                    continue                
                if b.pos < a.aend:
                    continue
                # todo: validate based on the amount of clipping/matches?
                deletion_intervals.append(pybedtools.Interval(samfile.getrname(a.tid), a.aend, b.pos, name=qname + "-" +  str(qname_idx)))
    return deletion_intervals

def validate_deletions(predicted_deletion, validated_deletions):
    return pybedtools.BedTool([predicted_deletion]).intersect(validated_deletions, wao=True)

def complement(s): 
    """Return the complementary sequence string.""" 
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    letters = list(s) 
    letters = [basecomplement[base] for base in letters] 
    return ''.join(letters) 
     
def reversecomplement(s): 
    """Return the reverse complement of the dna string.""" 
    s = complement(s) 
    return s[::-1]

def run_command(command):
    sys.stderr.write(command + "\n")
    subprocess.call(command, shell=True)

def sga_assemble(temp_dir):
    sys.stderr.write("sga preprocess -o {0}/preprocess.out -p 1 {0}/flank_reads_f1.fq {0}/flank_reads_fq2.fq\n".format(temp_dir))
    subprocess.call( "sga preprocess -o {0}/preprocess.out -p 1 {0}/flank_reads_f1.fq {0}/flank_reads_fq2.fq".format(temp_dir), shell=True)
    sys.stderr.write("sga index {0}/preprocess.out\n".format(temp_dir))
    subprocess.call( "sga index {0}/preprocess.out".format(temp_dir), shell=True)
    sys.stderr.write("sga correct {0}/preprocess.out -k 11 -x 2 -d 128 -o {0}/preprocess.ec.fa\n".format(temp_dir))
    subprocess.call( "sga correct {0}/preprocess.out -k 11 -x 2 -d 128 -o {0}/preprocess.ec.fa".format(temp_dir), shell=True)
    sys.stderr.write("sga index {0}/preprocess.ec.fa\n".format(temp_dir))
    subprocess.call( "sga index {0}/preprocess.ec.fa".format(temp_dir), shell=True)
    sys.stderr.write("sga rmdup {0}/preprocess.ec.fa\n".format(temp_dir))
    subprocess.call( "sga rmdup {0}/preprocess.ec.fa".format(temp_dir), shell=True)
    sys.stderr.write("sga overlap {0}/preprocess.ec.rmdup.fa -m 12\n".format(temp_dir))
    subprocess.call( "sga overlap {0}/preprocess.ec.rmdup.fa -m 12".format(temp_dir), shell=True)
    sys.stderr.write("sga assemble {0}/preprocess.ec.rmdup.asqg.gz -m 15 -d 0 -g 0 -b 0 -l 100\n".format(temp_dir))
    subprocess.call( "sga assemble {0}/preprocess.ec.rmdup.asqg.gz -m 15 -d 0 -g 0 -b 0 -l 100".format(temp_dir), shell=True)
    return "{}/default-contigs.fa".format(temp_dir)

def spades_assemble(temp_dir):
    run_command("/g/whelanch/software/SPAdes-2.5.0-Linux/bin/spades.py --careful --pe1-1 {0}/flank_reads_f1.fq --pe1-2 {0}/flank_reads_fq2.fq -t 1 -o {0}".format(temp_dir))
    return "{0}/contigs.fasta".format(temp_dir)

def assemble_contigs(prediction, bamfile, bwa_ref, unmapped_reads_db, flank = 1000, keep_temp_dirs = False):
    temp_dir = tempfile.mkdtemp()
    records_by_name = None
    flanks_only = True
    try:
        saved_path = os.getcwd()
        os.chdir(temp_dir)
        pred_left = prediction.start
        flank_left = pred_left - flank
        pred_right = prediction.end
        flank_right = pred_right + flank
        left_region = prediction.chrom + ":" + str(flank_left) + "-" + str(pred_left)
        right_region = prediction.chrom + ":" + str(pred_right) + "-" + str(flank_right)
        sys.stderr.write("samtools view -H {0} > {1}/header.sam\n".format(bamfile, temp_dir))
        subprocess.call("samtools view -H {0} > {1}/header.sam".format(bamfile, temp_dir), shell=True)

        if flanks_only:
            sys.stderr.write("samtools view {0} {1} > {2}/left_flank_mappings.sam\n".format(bamfile, left_region, temp_dir))
            subprocess.call("samtools view {0} {1} > {2}/left_flank_mappings.sam".format(bamfile, left_region, temp_dir), shell=True)
            sys.stderr.write("samtools view {0} {1} > {2}/right_flank_mappings.sam\n".format(bamfile, right_region, temp_dir))
            subprocess.call( "samtools view {0} {1} > {2}/right_flank_mappings.sam".format(bamfile, right_region, temp_dir), shell=True)
            run_command("cat {0}/left_flank_mappings.sam {0}/right_flank_mappings.sam | cut -f 1 | sort -u > {0}/flank_read_names.txt".format(temp_dir))
        else:
            full_region = prediction.chrom + ":" + str(flank_left) + "-" + str(flank_right)
            run_command("samtools view {0} {1} > {2}/region_mappings.sam\n".format(bamfile, full_region, temp_dir))
            run_command("cat {0}/region_mappings.sam | cut -f 1 | sort -u > {0}/flank_read_names.txt".format(temp_dir))

        unmapped_reads_file = open("{0}/unmapped_mates.sam".format(temp_dir), "w")
        c = unmapped_reads_db.cursor()
        for read_name in open("{0}/flank_read_names.txt".format(temp_dir)):
            sys.stderr.write("fetching results for read {0}\n".format(read_name.strip()))
            c.execute("SELECT sam FROM mappings WHERE readid = ?", [read_name.strip()])
            res=c.fetchone()
            if res is None:
                continue
            sys.stderr.write("got a read for read id {0}: {1}\n".format(read_name, res[0]))
            unmapped_reads_file.write(res[0].strip() + "\n")
        unmapped_reads_file.close()
        c.close()

#    run_command("samtools view -f 0x4 {0} | grep -F -f {1}/flank_read_names.txt > {1}/unmapped_mates.sam".format(bamfile, temp_dir))

        if flanks_only:
            run_command("cat {0}/left_flank_mappings.sam {0}/right_flank_mappings.sam {0}/unmapped_mates.sam | sort -k1,1 > {0}/mappings.sam".format(temp_dir))
        else:
            run_command("cat {0}/region_mappings.sam {0}/unmapped_mates.sam | sort -k1,1 > {0}/mappings.sam".format(temp_dir))
    
        reads_by_name = {}
        for read in open("{0}/mappings.sam".format(temp_dir), "r"):
            read_name = read.split("\t")[0]
            if not read_name in reads_by_name:
                reads_by_name[read_name] = []
            reads_by_name[read_name].append(read)
        filtered_mappings = open("{0}/mappings_filtered.sam".format(temp_dir), "w")
        for read_name in sorted(reads_by_name.keys()):
            if len(reads_by_name[read_name]) != 2:
                continue
            filtered_mappings.write(reads_by_name[read_name][0])
            filtered_mappings.write(reads_by_name[read_name][1])
        filtered_mappings.close()
        
        

    # sys.stderr.write("cat {0}/left_flank_read_names.txt {0}/right_flank_read_names.txt | sort -u > {0}/flank_read_names.txt\n".format(temp_dir))
    # subprocess.call( "cat {0}/left_flank_read_names.txt {0}/right_flank_read_names.txt | sort -u > {0}/flank_read_names.txt".format(temp_dir), shell=True)
    # sys.stderr.write("samtools view {0} | grep -F -f {1}/flank_read_names.txt | sort -k1,1 >> {1}/flank_reads.sam\n".format(bamfile, temp_dir))
    # subprocess.call("samtools view {0} | grep -F -f {1}/flank_read_names.txt | sort -k1,1 >> {1}/flank_reads.sam".format(bamfile, temp_dir), shell=True)
    
    # kmer_file = open("{0}/kmers.txt".format(temp_dir), "w")
    # kmer_size = 32
    # sys.stderr.write("looking for kmers of size {0}...\n".format(kmer_size))
    # for line in open("{0}/flank_reads.sam".format(temp_dir), "r"):
    #     if line.startswith("@"):
    #         continue
    #     seq = line.split("\t")[9]
    #     kmers = set()
    #     for i in range(0, len(seq) - kmer_size):
    #         kmer = seq[i:i+kmer_size]
    #         kmers.add(kmer)
    #         kmers.add(reversecomplement(kmer))
    #     kmer_file.write("\n".join(kmers))
    # kmer_file.close()
            
    # sys.stderr.write("samtools view -F 0x2 {0} | grep -F -f {1}/kmers.txt | cut -f1 | sort -u > {1}/reads_with_kmers.txt\n".format(bamfile, temp_dir))
    # subprocess.call("samtools view -F 0x2 {0} | grep -F -f {1}/kmers.txt | cut -f1 | sort -u > {1}/reads_with_kmers.txt\n".format(bamfile, temp_dir), shell=True)        

    # sys.stderr.write("samtools view {0} | grep -F -f {1}/reads_with_kmers.txt | sort -k1,1 >> {1}/flank_reads.sam\n".format(bamfile, temp_dir))
    # subprocess.call("samtools view {0} | grep -F -f {1}/reads_with_kmers.txt | sort -k1,1 >> {1}/flank_reads.sam\n".format(bamfile, temp_dir), shell=True)

        sys.stderr.write("cat {0}/header.sam {0}/mappings_filtered.sam | sort -u | samtools view -Sb - | samtools sort -n - {0}/flank_reads\n".format(temp_dir))   
        subprocess.call( "cat {0}/header.sam {0}/mappings_filtered.sam | sort -u | samtools view -Sb - | samtools sort -n - {0}/flank_reads".format(temp_dir), shell=True)
        sys.stderr.write("bamToFastq -i {0}/flank_reads.bam -fq {0}/flank_reads_f1.fq -fq2 {0}/flank_reads_fq2.fq\n".format(temp_dir))
        subprocess.call( "bamToFastq -i {0}/flank_reads.bam -fq {0}/flank_reads_f1.fq -fq2 {0}/flank_reads_fq2.fq".format(temp_dir), shell=True)

        contigs_file = spades_assemble(temp_dir)

        run_command("bwa mem {0} {1} > {2}/bwa_mem_out.sam\n".format(bwa_ref, contigs_file, temp_dir))
        (samfile, records_by_name) = read_sam_file("{0}/bwa_mem_out.sam".format(temp_dir))
    finally:
        os.chdir(saved_path)
        if not keep_temp_dirs:
            shutil.rmtree(temp_dir)

    return (samfile, records_by_name)
                                               

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('deletion_preds_file')
    parser.add_argument('--keep_unmapped_reads_db', default=False, dest='keep_unmapped_reads_db', action='store_true')
    parser.add_argument('--keep_temp_dirs', default=False, dest='keep_temp_dirs', action='store_true')
    parser.add_argument('--use_existing_unmapped_reads_db')
    args = parser.parse_args()

    print args 

    deletion_preds_file = args.deletion_preds_file
    keep_unmapped_reads_db = args.keep_unmapped_reads_db
    keep_temp_dirs = args.keep_temp_dirs

    bam_file = "/l2/users/whelanch/gene_rearrange/sv/jcvi_sim_chr2_allindels_100bp_dip/human_b36_male_chr2_venterindels_c15_i100_s30_rl100_sort.bam"
    unmapped_reads_db_name = tempfile.NamedTemporaryFile().name
    try:
        unmapped_reads_db = None
        if args.use_existing_unmapped_reads_db is None:
            sys.stderr.write("Creating unmapped reads database {0}..\n".format(unmapped_reads_db_name))
            start_create_db_time = int(round(time.time() * 1000))
            unmapped_reads_db = sqlite3.connect(unmapped_reads_db_name)
            c = unmapped_reads_db.cursor()
            c.execute('''CREATE TABLE mappings (readid text, sam text)''')
            view_process = subprocess.Popen("samtools view -f 0x4 {0}".format(bam_file), shell=True, stdout=subprocess.PIPE)
            i = 0
            time_cp = int(round(time.time() * 1000))
            for line in view_process.stdout:
                read_name = line[0:line.index("\t")]
                c.execute("insert into mappings values (?, ?)", (read_name, line.strip()))
                i = i + 1
                if i % 100000 == 0:
                    current_time = int(round(time.time() * 1000))
                    sys.stderr.write("loaded {0} records in {1}\n".format(i, (current_time - time_cp)))
                    time_cp = current_time
                    sys.stderr.write("wrote {0} reads to {1}\n".format(i, unmapped_reads_db_name))
                    unmapped_reads_db.commit()
            sys.stderr.write("creating index..")
            c.execute('''CREATE INDEX map_read_id_idx on mappings (readid)''')
            current_time = int(round(time.time() * 1000))
            sys.stderr.write("done indexing.. created db in {0}\n".format((current_time - start_create_db_time)))
        else:
            unmapped_reads_db_name = args.use_existing_unmapped_reads_db
            sys.stderr.write("using existing unmapped reads db {0}\n".format(unmapped_reads_db_name))
            unmapped_reads_db = sqlite3.connect(unmapped_reads_db_name)
        preds = pybedtools.BedTool(deletion_preds_file)
        results_file = open("results.txt", "w")
        for prediction in preds:
            sys.stderr.write("validating deletion pred " + str(prediction) + "\n")
            (samfile, records_by_name) = assemble_contigs(prediction, bam_file, "/l2/users/whelanch/genome_refs/10KG/hg18/human_b36_male_chr2.fasta", unmapped_reads_db, flank=1000, keep_temp_dirs=keep_temp_dirs)
            deletions_from_contigs = get_deletion_intervals(samfile, records_by_name)
            sys.stderr.write("\n".join(map(str,deletions_from_contigs)))
            sys.stderr.write("\n")
            if (len(deletions_from_contigs) == 0):
                sys.stderr.write("didn't get any deletions from the assembly\n")
                results_file.write(str(prediction) + "\tno assembled deletions\n")
                continue
            sys.stderr.write("comparing to truth\n")
            deletion_overlap = str(validate_deletions(prediction, deletions_from_contigs))
            sys.stderr.write("deletion_overlap: " +  deletion_overlap)
            results_file.write(deletion_overlap)

        results_file.close()
    finally:
        if args.use_existing_unmapped_reads_db is None:
            if keep_unmapped_reads_db:
                sys.stderr.write("Unmapped reads DB: {0}\n".format(unmapped_reads_db_name))
            else:
                sys.stderr.write("Removing temporary unmapped reads database..\n")
                os.remove(unmapped_reads_db_name)
    sys.stderr.write("Done.")

if __name__ == "__main__":
    main()
