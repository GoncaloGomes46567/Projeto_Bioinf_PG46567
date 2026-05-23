Este script automatiza a conversão e formatação de dados exportados do \*\*QIIME2\*\* nos 3 ficheiros de texto estruturados (`otu\_table.txt`, `taxonomy.txt`, `metadata.txt`). Após tratamento pelo script, os dados estão prontos para o upload na plataforma \*\*MicrobiomeAnalyst\*\*.









É necessário ter a biblioteca "pandas" instalada



Edita o ficheiro config/mapping.json para definir as regras do teu projeto (Mapeamento de Prefixo da Amostra -> Nome do Grupo).



Executa o script:



python scripts/prep\_data.py \\

&#x20; -i data/feature-table.tsv \\

&#x20; -f data/feature.fasta \\

&#x20; -t data/taxonomy.tsv \\

&#x20; -m config/mapping.json \\

&#x20; -o output



