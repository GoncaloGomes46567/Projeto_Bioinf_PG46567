import os
import argparse
import hashlib
import json
import pandas as pd

def carregar_mapa_fasta(caminho_fasta):
    mapa = {}
    with open(caminho_fasta, 'r') as f:
        id_atual, partes = None, []
        for linha in f:
            linha = linha.strip()
            if linha.startswith(">"):
                if id_atual and partes:
                    chave = hashlib.md5("".join(partes).encode()).hexdigest()
                    mapa[chave] = id_atual
                id_atual = linha[1:].split()[0]
                partes = []
            else:
                partes.append(linha)
        if id_atual and partes:
            chave = hashlib.md5("".join(partes).encode()).hexdigest()
            mapa[chave] = id_atual
    return mapa

def limpar_taxonomia(texto_tax):
    niveis   = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    prefixos = ["k__",     "p__",    "c__",   "o__",   "f__",    "g__",   "s__"]
    resultado = {nivel: "" for nivel in niveis}
    if pd.isna(texto_tax):
        return resultado
    for parte in str(texto_tax).split(";"):
        parte = parte.strip()
        for pref, nivel in zip(prefixos, niveis):
            if parte.startswith(pref):
                resultado[nivel] = parte[len(pref):].strip()
    return resultado

def descobrir_grupo(nome_amostra, dicionario_grupos):
    """Determina o grupo pelo prefixo mais longo (V.CN antes de V)."""
    for prefixo in sorted(dicionario_grupos.keys(), key=len, reverse=True):
        if nome_amostra.startswith(prefixo):
            return dicionario_grupos[prefixo]
    return "Unknown"

def main():
    parser = argparse.ArgumentParser(
        description="Prepara dados do QIIME2 para upload no MicrobiomeAnalyst."
    )
    parser.add_argument("-i", "--tabela-input", required=True,
                        help="Tabela de abundâncias TSV do QIIME2 (contagens brutas).")
    parser.add_argument("-f", "--fasta",        required=True,
                        help="Ficheiro FASTA com as sequências ASV (feature.fasta).")
    parser.add_argument("-t", "--taxonomia",    required=True,
                        help="Ficheiro de taxonomia TSV (feature.tax_assignments.txt).")
    parser.add_argument("-m", "--mapeamento",   required=True,
                        help='JSON com prefixo->grupo. Ex: {"V.CN":"NC","V":"Inoculum"}')
    parser.add_argument("-o", "--pasta-saida",  default="output_microbiome",
                        help="Pasta de destino para os ficheiros de output.")
    args = parser.parse_args()

    os.makedirs(args.pasta_saida, exist_ok=True)

    with open(args.mapeamento, 'r') as f:
        dados_grupos = json.load(f)

    print("-> Carregando tabela de abundâncias...")
    df_otu = pd.read_csv(args.tabela_input, sep="\t", skiprows=1, index_col=0)

    mapa_ids = carregar_mapa_fasta(args.fasta)
    df_otu.index = [mapa_ids.get(idx, idx) for idx in df_otu.index]

    df_otu = df_otu.astype(int)
    df_otu.index.name = "#NAME"
    df_otu.to_csv(os.path.join(args.pasta_saida, "otu_table.txt"), sep="\t")
    print(f"   otu_table.txt      ({df_otu.shape[0]} ASVs x {df_otu.shape[1]} amostras) — contagens brutas")

    otu_tss = df_otu.div(df_otu.sum(axis=0), axis=1)
    desvio_max = (otu_tss.sum(axis=0) - 1.0).abs().max()
    assert desvio_max < 1e-9, f"Validação TSS falhou — desvio máximo: {desvio_max}"
    otu_tss.index.name = "#NAME"
    otu_tss.to_csv(os.path.join(args.pasta_saida, "otu_table_tss.txt"), sep="\t", float_format="%.15g")
    print(f"   otu_table_tss.txt  ({df_otu.shape[0]} ASVs x {df_otu.shape[1]} amostras) — TSS normalizado, desvio máx={desvio_max:.2e}")

    print("-> Formatando taxonomia...")
    df_tax = pd.read_csv(args.taxonomia, sep="\t", index_col=0)
    df_tax.index = [mapa_ids.get(idx, idx) for idx in df_tax.index]

    tax_limpa = df_tax[df_tax.columns[0]].apply(limpar_taxonomia).apply(pd.Series)
    tax_limpa = tax_limpa.reindex(df_otu.index).fillna("")
    tax_limpa.index.name = "#TAXONOMY"
    tax_limpa.to_csv(os.path.join(args.pasta_saida, "taxonomy.txt"), sep="\t")
    print(f"   taxonomy.txt   ({tax_limpa.shape[0]} ASVs, {tax_limpa.shape[1]} níveis)")

    print("-> Gerando metadados...")
    lista_grupos = [descobrir_grupo(a, dados_grupos) for a in df_otu.columns]

    desconhecidos = [s for s, g in zip(df_otu.columns, lista_grupos) if g == "Unknown"]
    if desconhecidos:
        print(f"   AVISO: {len(desconhecidos)} amostra(s) sem grupo: {desconhecidos}")

    df_meta = pd.DataFrame({"Group": lista_grupos}, index=df_otu.columns)
    df_meta.index.name = "#NAME"
    df_meta.to_csv(os.path.join(args.pasta_saida, "metadata.txt"), sep="\t")
    print(f"   metadata.txt   ({df_meta.shape[0]} amostras, {df_meta['Group'].nunique()} grupos)")

    print(f"\nFicheiros guardados em: '{args.pasta_saida}/'")
    print("     MicrobiomeAnalyst: Taxonomy labels=QIIME | Normalized data=NAO marcar")
if __name__ == "__main__":
    main()
