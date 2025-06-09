
import argparse
from Bio import SeqIO

def filtrar_diamond(input_diamond, min_coverage=0.4, min_identity=40):
    cds_filtrados = set()
    with open(input_diamond, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            qseqid = fields[0]
            qlen = int(fields[1])
            length = int(fields[5])
            pident = float(fields[6])

            coverage = length / qlen
            if coverage >= min_coverage and pident >= min_identity:
                cds_filtrados.add(qseqid)
    return cds_filtrados

def extrair_sequencias_prodigal(input_prodigal, cds_filtrados):
    sequencias = {}
    for record in SeqIO.parse(input_prodigal, 'fasta'):
        if record.id in cds_filtrados:
            sequencias[record.id] = str(record.seq)
    return sequencias

def gerar_lista_filtrada(diamond_file, prodigal_file, output_file):
    cds_filtrados = filtrar_diamond(diamond_file)
    sequencias = extrair_sequencias_prodigal(prodigal_file, cds_filtrados)

    with open(output_file, 'w') as out:
        for cds_id, seq in sequencias.items():
            out.write(f">{cds_id}\n{seq}\n")

def main():
    parser = argparse.ArgumentParser(description="Filtrar CDS alinhados e extrair sequências do Prodigal")
    parser.add_argument("diamond_file", help="Arquivo de output do Diamond")
    parser.add_argument("prodigal_file", help="Arquivo de proteínas do Prodigal")
    parser.add_argument("output_file", help="Arquivo de saída contendo os CDS filtrados")

    args = parser.parse_args()
    gerar_lista_filtrada(args.diamond_file, args.prodigal_file, args.output_file)

if __name__ == "__main__":
    main()

