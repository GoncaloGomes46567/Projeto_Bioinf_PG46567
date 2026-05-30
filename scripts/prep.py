import os
import argparse
import hashlib
import json
import pandas as pd

def build_id_map(fasta_path):
    """Mapeia o hash MD5 da sequência para o ID da ASV."""
    mapping = {}
    with open(fasta_path, 'r') as f:
        name, parts = None, []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name and parts:
                    mapping[hashlib.md5("".join(parts).encode()).hexdigest()] = name
                name, parts = line[1:], []
            else:
                parts.append(line)
        if name and parts:
            mapping[hashlib.md5("".join(parts).encode()).hexdigest()] = name
    return mapping

def parse_taxonomy(tax_str):
    """Divide a string do QIIME2/Greengenes em colunas taxonómicas."""
    levels = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    prefixes = ["k__", "p__", "c__", "o__", "f__", "g__", "s__"]
    result = {l: "" for l in levels}
    for part in str(tax_str).split(";"):
        part = part.strip()
        for pref, level in zip(prefixes, levels):
            if part.startswith(pref):
                result[level] = part[len(pref):].strip()
    return result

def get_group(sample, groups_dict):
    """Determina o grupo com base no prefixo mais longo."""
    for prefix in sorted(groups_dict.keys(), key=len, reverse=True):
        if sample.startswith(prefix):
            return groups_dict[prefix]
    return "Unknown"

def main():
    parser = argparse.ArgumentParser(description="Prepara dados do QIIME2 para o MicrobiomeAnalyst.")
    parser.add_argument("-i", "--input-table", required=True, help="Caminho para a tabela de abundâncias (.tsv)")
    parser.add_argument("-f", "--fasta", required=True, help="Caminho para o ficheiro FASTA (.fasta)")
    parser.add_argument("-t", "--taxonomy", required=True, help="Caminho para o ficheiro de taxonomia (.tsv)")
    parser.add_argument("-m", "--mapping", required=True, help="Caminho para o JSON de mapeamento de grupos")
    parser.add_argument("-o", "--out-dir", default="output", help="Pasta de saída para os resultados")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    with open(args.mapping, 'r') as f:
        groups_dict = json.load(f)

    print("1. Processando Tabela de OTUs...")
    df = pd.read_csv(args.input_table, sep="\t", skiprows=1, index_col=0)
    id_map = build_id_map(args.fasta)
    df.index = df.index.map(lambda x: id_map.get(x, x))
    

    otu_float = df.astype(float)
    otu_tss = otu_float.div(otu_float.sum(axis=0), axis=1)
    otu_tss.index.name = "#NAME"
    otu_tss.to_csv(os.path.join(args.out_dir, "otu_table.txt"), sep="\t", float_format="%.6f")

    print("2. Processando Taxonomia...")
    tax_df = pd.read_csv(args.taxonomy, sep="\t", index_col=0)
    tax_col = tax_df.columns[0]  # Assume a primeira coluna como a string taxonómica
    parsed_tax = tax_df[tax_col].apply(parse_taxonomy).apply(pd.Series)
    parsed_tax = parsed_tax.loc[parsed_tax.index.isin(otu_tss.index)].reindex(otu_tss.index)
    parsed_tax.index.name = "#TAXONOMY"
    parsed_tax.to_csv(os.path.join(args.out_dir, "taxonomy.txt"), sep="\t")

    print("3. Processando Metadados...")
    meta = pd.DataFrame({"Group": [get_group(s, groups_dict) for s in otu_tss.columns]}, index=otu_tss.columns)
    meta.index.name = "#NAME"
    meta.to_csv(os.path.join(args.out_dir, "metadata.txt"), sep="\t")

    print(f"\n Concluído com sucesso! Ficheiros guardados na pasta '{args.out_dir}'")

if __name__ == "__main__":
    main()