import cyvcf2 as cyvcf
import pandas as pd

def walk_bcf(bcf):

    reader = cyvcf.VCF(bcf, gts012 = True)
    variants = []

    AC_THRESHOLD = 1

    filter_list = []
    chrom_list = []
    pos_list = []
    ref_list = []
    alt_list = []
    ac_list = []
    af_list = []
    depth_list = []
    abe_list = []

    for v in reader:
        AC = v.INFO.get("AC")
        AN = v.INFO.get("AN")
        FILTER = "PASS" if v.FILTER is None else v.FILTER
        
        if AC >= AC_THRESHOLD:
            filter_list.append(FILTER)
            chrom = v.CHROM
            pos = v.POS
            ref = v.REF
            alt = "".join(v.ALT)
            ac_list.append(AC)
            AF = AC / AN
            af_list.append(AF)
            abe_list.append(v.INFO.get("ABE"))
            depth_list.append(v.INFO.get("AVGDP"))
            variants.append(f"{chrom}-{pos}-{ref}-{alt}")

            chrom_list.append(chrom)
            pos_list.append(pos)
            ref_list.append(ref)
            alt_list.append(alt)

    variant_meta = pd.DataFrame({
            "CHROM" : chrom_list,
            "POS"   : pos_list,
            "REF"   : ref_list,
            "ALT"   : alt_list,
            "FILTER": filter_list,
            "AC"    : ac_list,
            "AF"    : af_list,
            "ABE"   : abe_list,
            "DEPTH" : depth_list,
            "ID"    : variants
    })

    return variant_meta
