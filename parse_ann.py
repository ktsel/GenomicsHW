import sys

def parse_ann(ann_field):
    """
    Парсит поле ANN и возвращает список аннотаций.
    """
    annotations = []
    for ann in ann_field.split(','):
        fields = ann.split('|')
        annotation = {
            'Allele': fields[0],
            'Annotation': fields[1],
            'Annotation_Impact': fields[2],
            'Gene_Name': fields[3],
            'Gene_ID': fields[4],
            'Feature_Type': fields[5],
            'Feature_ID': fields[6],
            'Transcript_BioType': fields[7],
            'Rank': fields[8],
            'HGVS.c': fields[9],
            'HGVS.p': fields[10],
            'cDNA.pos / cDNA.length': fields[11],
            'CDS.pos / CDS.length': fields[12],
            'AA.pos / AA.length': fields[13],
            'Distance': fields[14],
            'ERRORS / WARNINGS / INFO': fields[15]
        }
        annotations.append(annotation)
    return annotations

def main(input_file, output_file):
    """
    Читает файл с аннотациями и записывает их в читаемом формате.
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Заголовок таблицы
        outfile.write("Allele\tAnnotation\tAnnotation_Impact\tGene_Name\tGene_ID\tFeature_Type\tFeature_ID\tTranscript_BioType\tRank\tHGVS.c\tHGVS.p\tcDNA.pos / cDNA.length\tCDS.pos / CDS.length\tAA.pos / AA.length\tDistance\tERRORS / WARNINGS / INFO\n")
        
        for line in infile:
            if line.startswith("ANN"):
                continue  # Пропускаем заголовок
            ann_field = line.strip()
            annotations = parse_ann(ann_field)
            for ann in annotations:
                outfile.write("\t".join(ann.values()) + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python parse_ann.py <input_file> <output_file>")
    else:
        main(sys.argv[1], sys.argv[2])