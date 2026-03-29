# 0. AUTO-INSTALLER
import os
try:
    import biopython
except ImportError:
    os.system('pip install biopython pandas')

from Bio import Entrez
import pandas as pd
import time

Entrez.email = "leafcloverfive@gmail.com"
Entrez.tool = "ISB_GWAS_PanAncestry_Miner"

print("\n" + "=" * 70)
print("INITIATING PAN-ANCESTRY GWAS DATA MINER (MDD: EAS, AFR, SAS, AMR)")
print("=" * 70)

query = '''
("Major Depressive Disorder"[Title/Abstract] OR MDD[Title/Abstract]) 
AND ("Genome-Wide Association"[Title/Abstract] OR GWAS[Title/Abstract] OR "Polygenic Risk"[Title/Abstract] OR PRS[Title/Abstract]) 
AND ("East Asian"[Title/Abstract] OR "African"[Title/Abstract] OR "South Asian"[Title/Abstract] OR "Hispanic"[Title/Abstract] OR EAS[Title/Abstract] OR AFR[Title/Abstract] OR SAS[Title/Abstract] OR AMR[Title/Abstract])
AND open access[filter]
'''

try:
    search_handle = Entrez.esearch(db="pmc", term=query, retmax=50)
    search_record = Entrez.read(search_handle)
    search_handle.close()

    id_list = search_record.get("IdList", [])
    
    if not id_list:
        print("Pencarian belum menemukan hasil.")
    else:
        results = []
        for pmcid in id_list:
            try:
                time.sleep(0.5) 
                summary_handle = Entrez.esummary(db="pmc", id=pmcid)
                summary_record = Entrez.read(summary_handle)
                summary_handle.close()
                
                record = summary_record[0] if isinstance(summary_record, list) and len(summary_record) > 0 else {}
                title = record.get("Title", "No Title Available")
                pub_date = record.get("PubDate", "Unknown Date")
                
                results.append({
                    "PMCID": f"PMC{pmcid}",
                    "Publication_Date": pub_date,
                    "Title": title
                })
                print(f"Data Genetik Terkunci: PMC{pmcid}")
            except Exception as inner_e:
                continue

        df_gwas = pd.DataFrame(results)
        df_gwas.to_csv('NCBI_GWAS_PanAncestry_MDD.csv', index=False)
        print("EKSTRAKSI DATA GENETIK GLOBAL SELESAI.")

except Exception as e:
    print(f"CRITICAL ERROR: {e}")
