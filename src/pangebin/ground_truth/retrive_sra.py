from __future__ import annotations

import time

from Bio import Entrez

Entrez.email = "mttsgro@gmail.com"
Entrez.tool = "plasbin-flow@pangenome"
Entrez.api_key = "7c5ce3b0327c240c68dc660082564a0bda08"
WAIT_TIME = 1


def sra_from_biosample(biosampleID: str):
    def src(database: str, term: str, retmax: int | None = 10) -> dict:
        """Searches entrez for a term and returns the results."""  # noqa: D401
        time.sleep(WAIT_TIME)
        handle = Entrez.esearch(db=database, retmax=retmax, term=term)
        record = Entrez.read(handle, validate=False)
        handle.close()
        return record

    def summ(database: str, id: str, retmax: int | None = 10) -> dict:
        """Retrieves a summary of an assembly from entrez."""  # noqa: D401
        time.sleep(WAIT_TIME)
        handle = Entrez.esummary(db=database, retmax=retmax, id=id)
        record = Entrez.read(handle, validate=False)
        handle.close()
        return record

    biosamples = src("BioSample", biosampleID)
    for x in biosamples["IdList"]:
        bio_summ = summ("BioSample", x)
        sample_data = bio_summ["DocumentSummarySet"]["DocumentSummary"][0]["SampleData"]
        str_start = sample_data.find('<Id db="SRA">')
        str_end = sample_data.find("</Ids>")
        offset = 1 + sample_data[str_start:str_end].find(">")
        str_start += offset
        str_end = str_start + sample_data[str_start:str_end].find("<")
        srs = sample_data[str_start:str_end]
        sra_ids = src("SRA", srs)
        for y in sra_ids["IdList"]:
            sra_summ = summ("SRA", y)
            sra_info = str(sra_summ[0]["ExpXml"])
            if str(sra_info).lower().find("illumina") > 0:
                print(sra_summ)
                print(y, "its ILLUMINA!")
                print("SRS:", srs)
                srx_off = sra_info.find("SRX")
                srx_len = sra_info[sra_info.find("SRX") :].find('"')
                srx = sra_info[srx_off : srx_off + srx_len]
                print("SRX:", srx)
                sra_off = sra_info.find("SRA")
                sra_len = sra_info[sra_info.find("SRA") :].find('"')
                sra = sra_info[sra_off : sra_off + sra_len]
                print("SRA:", sra)
            else:
                print(y, "no illumina, RIP")


def main():
    sra_from_biosample("SAMN16357463")
