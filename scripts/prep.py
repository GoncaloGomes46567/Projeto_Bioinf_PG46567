import os
import argparse
import hashlib
import json
import pandas as pd

def carregar_mapa_fasta(caminho_fasta):
    mapa = {}
    with open(caminho_fasta, 'r') as f:
        id_atual = None
        for linha in f:
            linha = linha.strip()
            if linha.startswith(">"):
                id_atual = linha[1:].split()[0]
                mapa[id_atual] = id_atual
    return mapa
    
def limpar_taxonomia(texto_tax):
    niveis = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    prefixos = ["k__", "p__", "c__", "o__", "f__", "g__", "s__"]
    resultado = {nivel: "" for nivel in niveis}
    if pd.isna(texto_tax):
        return resultado
    partes = str(texto_tax).split(";")
    for parte in partes:
        parte = parte.strip()
        for pref, nivel in zip(prefixos, niveis):
            if parte.startswith(pref):
                resultado[nivel] = parte[len(pref):].strip()
    return resultado

def descobrir_grupo(nome_amostra, dicionario_grupos):
    prefixos_ordenados = sorted(dicionario_grupos.keys(), key=len, reverse=True)
    for prefixo in prefixos_ordenados:
        if nome_amostra.startswith(prefixo):
            return dicionario_grupos[prefixo]
    return "Unknown"

def main():
    parser = argparse.ArgumentParser(description="Prepara dados do QIIME2 formatados para o MicrobiomeAnalyst.")
    parser.add_argument("-i", "--tabela-input", required=True)
    parser.add_argument("-f", "--fasta", required=True)
    parser.add_argument("-t", "--taxonomia", required=True)
    parser.add_argument("-m", "--mapeamento", required=True)
    parser.add_argument("-o", "--pasta-saida", default="output_microbiome")
    args = parser.parse_args()

    os.makedirs(args.pasta_saida, exist_ok=True)

    with open(args.mapeamento, 'r') as f:
        dados_grupos = json.load(f)

    df_otu = pd.read_csv(args.tabela_input, sep="\t", skiprows=1, index_col=0)
    mapa_ids = carregar_mapa_fasta(args.fasta)
    df_otu.index = [mapa_ids.get(idx, idx) for idx in df_otu.index]
    
    otu_float = df_otu.astype(float)
    otu_tss = otu_float.div(otu_float.sum(axis=0), axis=1)
    
    otu_tss.index.name = "#NAME"
    caminho_otu = os.path.join(args.pasta_saida, "otu_table.txt")
    otu_tss.to_csv(caminho_otu, sep="\t", float_format="%.6f")

    df_tax = pd.read_csv(args.taxonomia, sep="\t", index_col=0)
    df_tax.index = [mapa_ids.get(idx, idx) for idx in df_tax.index]
    coluna_tax = df_tax.columns[0]
    tax_limpa = df_tax[coluna_tax].apply(limpar_taxonomia).apply(pd.Series)
    
    tax_limpa = tax_limpa.reindex(otu_tss.index)
    tax_limpa.index.name = "#TAXONOMY"
    caminho_tax = os.path.join(args.pasta_saida, "taxonomy.txt")
    tax_limpa.to_csv(caminho_tax, sep="\t")

    lista_grupos = []
    for amostra in otu_tss.columns:
        grupo = descobrir_grupo(amostra, dados_grupos)
        lista_grupos.append(grupo)
        
    df_meta = pd.DataFrame({"Group": lista_grupos}, index=otu_tss.columns)
    df_meta.index.name = "#NAME"
    caminho_meta = os.path.join(args.pasta_saida, "metadata.txt")
    df_meta.to_csv(caminho_meta, sep="\t")

if __name__ == "__main__":
    main()
