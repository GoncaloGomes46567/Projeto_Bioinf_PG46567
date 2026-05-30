import os
import argparse
import hashlib
import json
import pandas as pd

def carregar_mapa_fasta(caminho_fasta):
    """
    Lê o ficheiro FASTA e cria um mapa entre o ID da sequência 
    e o cabeçalho completo.
    """
    mapa = {}
    with open(caminho_fasta, 'r') as f:
        id_atual = None
        for linha in f:
            linha = linha.strip()
            if linha.startswith(">"):
                # Remove o '>' e guarda o ID (geralmente a primeira palavra)
                id_atual = linha[1:].split()[0]
                mapa[id_atual] = id_atual
    return mapa
    
def limpar_taxonomia(texto_tax):
    """
    Divide a string de taxonomia do QIIME2 (k__...; p__...)
    nas colunas certas que o MicrobiomeAnalyst pede.
    """
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
    """
    Descobre o grupo da amostra com base no prefixo mais longo 
   correspondente no ficheiro JSON.
    """
    prefixos_ordenados = sorted(dicionario_grupos.keys(), key=len, reverse=True)
    for prefixo in prefixos_ordenados:
        if nome_amostra.startswith(prefixo):
            return dicionario_grupos[prefixo]
    return "Unknown"

def main():
    parser = argparse.ArgumentParser(description="Prepara dados do QIIME2 formatados para o MicrobiomeAnalyst.")
    parser.add_argument("-i", "--tabela-input", required=True, help="Tabela de abundâncias (.tsv)")
    parser.add_argument("-f", "--fasta", required=True, help="Ficheiro FASTA (.fasta)")
    parser.add_argument("-t", "--taxonomia", required=True, help="Ficheiro de taxonomia (.tsv)")
    parser.add_argument("-m", "--mapeamento", required=True, help="JSON com o mapeamento dos grupos/metadados")
    parser.add_argument("-o", "--pasta-saida", default="output_microbiome", help="Pasta para guardar os resultados")
    args = parser.parse_args()

    os.makedirs(args.pasta_saida, exist_ok=True)


    with open(args.mapeamento, 'r') as f:
        dados_grupos = json.load(f)

    print("-> Formatando a tabela de abundâncias...")
    
<<<<<<< HEAD

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
=======
    df_otu = pd.read_csv(args.tabela_input, sep="\t", skiprows=1, index_col=0)
    mapa_ids = carregar_mapa_fasta(args.fasta)
    df_otu.index = [mapa_ids.get(idx, idx) for idx in df_otu.index]
    df_otu = df_otu.astype(int)
    df_otu.index.name = "#NAME"
    caminho_otu = os.path.join(args.pasta_saida, "otu_table.txt")
    df_otu.to_csv(caminho_otu, sep="\t")
    print("-> Formatando os dados taxonómicos...")
    df_tax = pd.read_csv(args.taxonomia, sep="\t", index_col=0)
    df_tax.index = [mapa_ids.get(idx, idx) for idx in df_tax.index]
    coluna_tax = df_tax.columns[0]
    tax_limpa = df_tax[coluna_tax].apply(limpar_taxonomia).apply(pd.Series)
    tax_limpa = tax_limpa.reindex(df_otu.index)
    tax_limpa.index.name = "#TAXONOMY"
    caminho_tax = os.path.join(args.pasta_saida, "taxonomy.txt")
    tax_limpa.to_csv(caminho_tax, sep="\t")
    print("-> Gerando ficheiro de metadados...")
    lista_grupos = []
    
    for amostra in df_otu.columns:
        grupo = descobrir_grupo(amostra, dados_grupos)
        lista_grupos.append(grupo)
    df_meta = pd.DataFrame({"Group": lista_grupos}, index=df_otu.columns)
    df_meta.index.name = "#NAME"
    caminho_meta = os.path.join(args.pasta_saida, "metadata.txt")
    df_meta.to_csv(caminho_meta, sep="\t")
    print(f"\n[OK] Tudo pronto! Os ficheiros foram guardados em: '{args.pasta_saida}'")
>>>>>>> 70c8f5725d8e4171dce7b1f32704bb63bef8ad08

if __name__ == "__main__":
    main()
