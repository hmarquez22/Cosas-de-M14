import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Entrez

Entrez.email = 'hmarquez22@ilg.cat'

def crear_record(cadena):
   handle = Entrez.esearch(db="nucleotide", term=cadena, retmax=2)
   recordintern = Entrez.read(handle)
   with open("brca2.fasta", "w") as resultats:
       for i in recordintern["IdList"]:
           handle1 = Entrez.efetch(db="nucleotide", id=i, rettype="fasta")
           resultats.write(handle1.read())

generic = "brca2 transcript variant"
record = crear_record(generic)

for i in SeqIO.parse("brca2.fasta","fasta"):
   handle = Entrez.efetch(db="nucleotide", id=i, rettype="gb")

   cont = 0
   for req in SeqIO.parse(handle,"genbank"):
       for feature in req.features:
           if feature == "exon":
               cont = cont +1
           if feature == "CDS":
               feature.location.extract(req).seq[150]
           if feature == "gene":

           if feature == "misc_feature"

   #valor = trencaid(i.id) #en valor vull només una part del id
   #diccionari[valor] = i.seq[:200]

     ##  claus = list(diccionario.keys())
  ## for i in range(len(claus)):
    ##   for j in range(i + 1, len(claus)):
         ##  simple_handle.write("-----------------------------------------------\n")
       ##    simple_handle.write("alineaments entre: " + claus[i] + "i" + claus[j] + "\n")
         ##  simple_handle.write("-----------------------------------------------\n")
       ##    valor = alinea(diccionario[claus[i]], diccionario[claus[j]])
        ##   simple_handle.write(str(diccionario[claus[i]]) + "\n")
      ##    simple_handle.write(str(diccionario[claus[j]]) + "\n")
         ##  simple_handle.write(str(valor) + "\n")
