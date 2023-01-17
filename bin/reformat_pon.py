#!/usr/bin/env python3

from pysam import VariantFile
import os

def reformat_pon(pon_file, out):
    with VariantFile(pon_file) as fr:
        header = fr.header
        header.info.add("PANEL_OF_NORMALS", number=0, type="Flag", description="Overlap with germline call among panel of normals")
        with VariantFile("tmp_.vcf", "w", header=header) as fw:
            for record in fr:
                record.info["set"] = "PANEL_OF_NORMALS"
                fw.write(record)

    os.system(f"mv tmp_.vcf {out}")
    os.system(f"bgzip {out}")
    os.system(f"tabix {out}.gz")

if __name__ == "__main__":
    from xcmds import xcmds

    xcmds.xcmds(locals())
