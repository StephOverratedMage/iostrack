# scripts/update_taxa_list.py
import io
import pandas as pd
import requests

def update_viral_taxa():
    # Queries ENA for all viral species (tax_tree(10239) = Viruses)
    url = "https://www.ebi.ac.uk/ena/portal/api/search"
    params = {
        "result": "taxon",
        "query": "tax_tree(10239) AND rank=\"species\"",
        "fields": "tax_id,scientific_name",
        "limit": 0,  # 0 = fetch all matching species
        "format": "tsv"
    }
    print("Fetching viral species list from ENA...")
    res = requests.get(url, params=params, timeout=60)
    res.raise_for_status()
    
    df = pd.read_csv(io.StringIO(res.text), sep="\t")
    df = df.rename(columns={"scientific_name": "Organism", "tax_id": "Taxon_ID"})
    df = df.sort_values("Organism").drop_duplicates(subset=["Organism"])
    
    df.to_csv("config/pathogens.csv", index=False)
    print(f"✅ Saved {len(df)} viral species to config/pathogens.csv!")

if __name__ == "__main__":
    update_viral_taxa()