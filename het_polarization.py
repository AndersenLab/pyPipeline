'''
Heterozygote Polarization Script
usage:
bcftools view <filename> | python het_polarization.py  | bcftools view -O b > <filename.het.polarized.bcf>

Tags variants 'pushed' to ref or alt as follows:

AA - Pushed towards reference
AB - Kept as het
BB - Pushed towards alternative
'''
 
import sys

het_prior = 0.000001
 
def phred2p(phred):
    return 10**(phred/-10.0)
 
def main():
    polarize_alt = 0
    polarize_ref = 0
    remain_het = 0
    format_added = False
    for l in sys.stdin.xreadlines():
        if l.startswith("#CHROM"):
            # Get Sample information and count
            samples = l.replace("\n","").split("\t")[9:]
            sys.stdout.write(l)
        elif l.startswith("#"):
            # Add Info line for het polarization flag
            if l.startswith("##FORMAT") and format_added == False:
                format_added = True
                sys.stdout.write("##FORMAT=<ID=HP,Number=1,Type=String,Description=\"Flag used to mark whether a variant was polarized\">\n")
            # Pass comment lines.
            sys.stdout.write(l)
        # Ignore indels, and see if there are any hets in the line before proceeding.
        elif l.find("0/1") != -1 and l.find("INDEL") == -1:
        # Check and see if there are any hets
            l = l.strip().split("\t")
            PL = l[8].split(":").index("PL")
            add_HP_flag = 0
            for k,v in enumerate(l[9:]):
                # Exclude any line 
                if v.startswith("0/1"):
                    PL_set = [phred2p(int(i)) for i in v.split(":")[PL].split(",")]
                    cond_prob = [PL_set[0]*((1-het_prior)/2), PL_set[1]*het_prior, PL_set[2]*(1-het_prior)/2]
                    cond_prob = [cond_prob[0]/sum(cond_prob), cond_prob[1]/sum(cond_prob), cond_prob[2]/sum(cond_prob)]
                    if add_HP_flag == 0:
                        l[8] = l[8] + ":HP"
                        add_HP_flag = 1
                    if (max(cond_prob) == cond_prob[0]):
                        l[k+9] = v.replace("0/1","0/0") + ":AA"   
                    elif (max(cond_prob) == cond_prob[2]):
                        l[k+9] = v.replace("0/1","1/1") + ":BB"
                    else:
                        l[k+9] = v + ":AB"
            sys.stdout.write("\t".join(l) + '\n')
        else:
            sys.stdout.write(l)

if __name__ == '__main__':
    main()


