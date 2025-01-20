from biopandas.pdb import PandasPdb
import re
import pymongo
from Bio import PDB
import Bio
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio import Align
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import pandas as pd
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import Phylo
from pylab import *
from Bio import motifs
import os
import pandas as pd
import matplotlib.pyplot as plt
import re
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from io import BytesIO
import base64
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Align
from Bio.pairwise2 import format_alignment


def biotraToolFunction(sequence):
        sequence=sequence.replace('*','')
        def isDNA(seq, alphabet='dna'):
             alphabets = {'dna': re.compile('^[acgtn]*$', re.I), 'protein': re.compile('^[acdefghiklmnpqrstvwy]*$', re.I)}
             if alphabets[alphabet].search(seq) is not None:
                return True
             else:
                return False
        
        results=[]
        sequence=sequence.upper()
        my_seq=Seq(sequence)
        if isDNA(sequence)==True:
             GC_quantity=gc_fraction(my_seq)#CANTIDAD DE GC  (Influye en la evolución de las proteínas por su coste energético, a mayor cantidad de GC menor costo para la síntesis de proteínas pero mayor costo para la síntesis de nucleótidos)
             A_quantity=sequence.count("A")/len(sequence)
             G_quantity=sequence.count("G")/len(sequence)
             C_quantity=sequence.count("C")/len(sequence)
             T_quantity=sequence.count("T")/len(sequence)
             GC_quantity=float("{:.4f}".format(GC_quantity))
             results.append(GC_quantity)#[0]
             complement=my_seq.complement() #COMPLEMENTARIA
             results.append(complement)#[1]
             reverse_Complement=my_seq.reverse_complement()#REVERSA COMPLEMENTARIA
             results.append(reverse_Complement)#[2]
             results.append(A_quantity)#[3]
             results.append(G_quantity)#[4]
             results.append(C_quantity)#[5]
             results.append(T_quantity)#[6]
             return(results)
        else:
             results=[]
             results.append("The sequence is not DNA")
             results.append("The sequence is not DNA")
             results.append("The sequence is not DNA")
             results.append("The sequence is not DNA")
             results.append("The sequence is not DNA")
             results.append("The sequence is not DNA")
             results.append("The sequence is not DNA")
             return(results)

def proteinAnalysisF(sequenceProts,pHValue):
     sequenceProts=sequenceProts.replace('*','')
     pHValue=float(pHValue)
     def isDNA(seq, alphabet='dna'):
             alphabets = {'dna': re.compile('^[acgtn]*$', re.I), 'protein': re.compile('^[acdefghiklmnpqrstvwy]*$', re.I)}
             if alphabets[alphabet].search(seq) is not None:
                  return True
             else:
                return False
             
     resultsProt=[]
     sequenceProts=sequenceProts.upper()
     sequenceProts=sequenceProts.replace("*","")
     if isDNA(sequenceProts)==False:
          sequenceProts=Seq(sequenceProts)
          analysis_Proteins=ProteinAnalysis(sequenceProts)
          resultsProt.append(analysis_Proteins)#0
          cantidad_aminoacidos=analysis_Proteins.count_amino_acids()
          
          resultsProt.append(cantidad_aminoacidos)#1
          cantidad_porcentual_aminoacidos=analysis_Proteins.get_amino_acids_percent()
          
          resultsProt.append(cantidad_porcentual_aminoacidos)#2
          peso_molecular=analysis_Proteins.molecular_weight()
          resultsProt.append(peso_molecular)#3
          aromaticity=analysis_Proteins.aromaticity()
          resultsProt.append(aromaticity)#4
          indice_inestabilidad=analysis_Proteins.instability_index()
          resultsProt.append(indice_inestabilidad)#5
          punto_isoelectrico=analysis_Proteins.isoelectric_point()
          resultsProt.append(punto_isoelectrico)#6
          estructura_secundaria=analysis_Proteins.secondary_structure_fraction()
          resultsProt.append(estructura_secundaria)#7
          epsilon_prot=analysis_Proteins.molar_extinction_coefficient()#reducidos y oxidados
          resultsProt.append(epsilon_prot)#8
          reducidos_cysteines=epsilon_prot[0]
          resultsProt.append(reducidos_cysteines)#9
          puentes_disulfid=epsilon_prot[1]
          resultsProt.append(puentes_disulfid)#10
          gravy=analysis_Proteins.gravy()
          resultsProt.append(gravy)#11
          #protein_scale=analysis_Proteins.protein_scale(window=9,param_dict=0)
          protein_scale=0
          resultsProt.append(protein_scale)#12
          flexibility=analysis_Proteins.flexibility()
          resultsProt.append(flexibility)#13
          carga_PH=analysis_Proteins.charge_at_pH(pHValue)
          resultsProt.append(carga_PH)#14
          resultsProt.append(len(sequenceProts))#15

          #GRAFICO Hidrofobicidad
          kd = { 'A': 1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C': 2.5,#kyte dolittle
                'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I': 4.5,
                'L': 3.8,'K':-3.9,'M': 1.9,'F': 2.8,'P':-1.6,
                'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V': 4.2 }
          values = []
          num_residues=len(sequenceProts)
          for residue in sequenceProts:
               values.append(kd[residue])
          x_data = range(1, num_residues+1)
          #ax = fig.add_axes([1, num_residues,0.75,0.75]) # axis starts at 0.1, 0.1
          #ax.set_title("K&D Hydrophobicity")
          #ax.set_xlabel("xlabel")
          #ax.plot(x_data, values)
          plt.plot(x_data,values)
          plt.xlabel="Residue Number"
          plt.ylabel="Hydrophobicity"
          
          
          xmin, xmax = xlim((1, num_residues))

          buf = BytesIO()
          
          plt.savefig(buf, format='jpg')
          image_base64 = base64.b64encode(buf.getvalue()).decode('utf-8').replace('\n', '')
          buf.close()
          resultsProt.append(image_base64)#16
          resultsProt.append(list(cantidad_aminoacidos.values()))#17
          resultsProt.append(cantidad_aminoacidos.keys())#18
          
          return(resultsProt)
     else:
          resultsProt=[]
          sequenceProts="ACGT"
          for i in range(13):
               resultsProt.append("This is not a Proteome")

          return(resultsProt)
          




def transcriptionFunction(DNAsequence):
     DNAsequence=DNAsequence.replace('*','')
     results=[]
     results.append(DNAsequence)#0-DNA sequence
     def isDNA(seq, alphabet='dna'):
          alphabets = {'dna': re.compile('^[acgtn]*$', re.I), 'protein': re.compile('^[acdefghiklmnpqrstvwy]*$', re.I)}
          if alphabets[alphabet].search(seq) is not None:
               return True
          else:
               return False
     if isDNA(DNAsequence)==True:
          codingDNA=Seq(DNAsequence)
          templateDNA=codingDNA.reverse_complement()
          results.append(templateDNA)#1-Template DNA
          transcriptionResult=codingDNA.transcribe()
          results.append(transcriptionResult)#2-Transcription Result of CODING DNA
          reverseDNATemplateStrand=templateDNA.reverse_complement()
          results.append(reverseDNATemplateStrand)#3-Reverse DNA Template Strand
          return(results)
     else:
          results.append("This is not a DNA sequence.")
          results.append("This is not a DNA sequence.")
          results.append("This is not a DNA sequence.")
          results.append("This is not a DNA sequence.")
          return(results)
     
def backTranscriptionTool(mRNASequence):
     def ismRNA(seq, alphabet='mRNA'):
          alphabets = {'dna': re.compile('^[acgtn]*$', re.I), 'protein': re.compile('^[acdefghiklmnpqrstvwy]*$', re.I),'mRNA': re.compile('^[acgu]*$', re.I)}
          if alphabets[alphabet].search(seq) is not None:
               return True
          else:
               return False
     results=[]
     if(ismRNA(mRNASequence))==True:
          results.append(mRNASequence)#0-mRNA original sequence
          mRNASequence=Seq(mRNASequence)
          backTranscriptionSequence=mRNASequence.back_transcribe()
          results.append(backTranscriptionSequence)#1-BackTranscriptionResult
          return(results)
     else:
          results.append("This is not a mRNA sequence.")
          results.append("This is not a mRNA sequence.")
          return(results)
     
##################################################################################################################################################

def translationTool(sequence,codonTable):
     sequence=sequence.replace('*','')
     results=[]
     results.append(sequence)#0-original sequence
     def ismRNA(seq, alphabet='mRNA'):
          alphabets = {'dna': re.compile('^[acgtn]*$', re.I), 'protein': re.compile('^[acdefghiklmnpqrstvwy]*$', re.I),'mRNA': re.compile('^[acgu]*$', re.I)}
          if alphabets[alphabet].search(seq) is not None:
               return True
          else:
               return False
          
     def isDNA(seq, alphabet='dna'):
          alphabets = {'dna': re.compile('^[acgtn]*$', re.I), 'protein': re.compile('^[acdefghiklmnpqrstvwy]*$', re.I)}
          if alphabets[alphabet].search(seq) is not None:
               return True
          else:
               return False
          
     def isProt(seq, alphabet='protein'):
          alphabets = {'dna': re.compile('^[acgtn]*$', re.I), 'protein': re.compile('^[acdefghiklmnpqrstvwy]*$', re.I)}
          if alphabets[alphabet].search(seq) is not None:
               return True
          else:
               return False
     if ismRNA(sequence)==True:
          results.append("This is a mRNA sequence")#1 type of sequence-mRNA
          stopCodon=""###############
          positions=3################
          while(positions > 0):######
               stopCodon=stopCodon+sequence[-positions]
               positions = positions-1
          mRNAsequence=Seq(sequence)
          if codonTable==100:
               listTranslations=[]
               listCodonTables=[1,2,3,4,5,6,9,10,11,12,13,14,15,16,21,22,23,24,25,26,27,28,29,30,31,33]
               for elemento in listCodonTables:
                    if elemento==11:
                         translationResult1=mRNAsequence.translate(table=elemento)
                         if len(sequence)%3==0 & (stopCodon=='TAA'or stopCodon=='TAG' or stopCodon=='TGA' ):
                              try:
                                   translationResult2=mRNAsequence.translate(table=11,cds=True)
                              except:
                                   translationResult2="Error"
                              cadena="If is not a CDS: "+translationResult1+" . If is a CDS: " + translationResult2
                         else:
                              translationResult2="This is not a valid Stop Codon"
                              cadena="If is not a CDS: "+translationResult1+" . If is a CDS: " + translationResult2


                         listTranslations.append(cadena)
                    else:
                         translationResult=mRNAsequence.translate(table=elemento)
                         listTranslations.append(translationResult)
               results.append(listTranslations)#2-lista de listas
          elif codonTable==11:
               translationResult1=mRNAsequence.translate(table=codonTable)
               if len(sequence)%3==0 & (stopCodon=='TAA'or stopCodon=='TAG' or stopCodon=='TGA' ):
                    try:
                         translationResult2=mRNAsequence.translate(table=11,cds=True)
                    except:
                         translationResult2="Error"
                    cadena="If is not a CDS: "+translationResult1+" . If is a CDS: " + translationResult2
               else:
                    translationResult2="This is not a valid Stop Codon"
                    cadena="If is not a CDS: "+translationResult1+" . If is a CDS: " + translationResult2
               
               results.append(cadena)#2 Translation Result
          else:
               translationResult=mRNAsequence.translate(table=codonTable)
               results.append(translationResult)#2 Translation Result
          return(results)


     elif isDNA(sequence)==True:
          results.append("This is a DNA sequence")#1 type of sequence-DNA
          stopCodon=""###############
          positions=3################
          while(positions > 0):######
               stopCodon=stopCodon+sequence[-positions]
               positions = positions-1
          DNAsequence=Seq(sequence)
          if codonTable==100:
               listTranslations=[]
               listCodonTables=[1,2,3,4,5,6,9,10,11,12,13,14,15,16,21,22,23,24,25,26,27,28,29,30,31,33]
               for elemento in listCodonTables:
                    if elemento==11:
                         translationResult1=DNAsequence.translate(table=elemento)
                         if len(sequence)%3==0 & (stopCodon=='TAA'or stopCodon=='TAG' or stopCodon=='TGA' ):
                              try:
                                   translationResult2=mRNAsequence.translate(table=11,cds=True)
                              except:
                                   translationResult2="Error"
                              cadena="If is not a CDS: "+translationResult1+" . If is a CDS: " + translationResult2
                         else:
                              translationResult2="This is not a valid Stop Codon"
                              cadena="If is not a CDS: "+translationResult1+" . If is a CDS: " + translationResult2
               
                         
                         listTranslations.append(cadena)
                    else:
                         translationResult=DNAsequence.translate(table=elemento)
                         listTranslations.append(translationResult)
               results.append(listTranslations)#2-lista de listas
          elif codonTable==11:
               translationResult1=DNAsequence.translate(table=codonTable)
               if len(sequence)%3==0 & (stopCodon=='TAA' or stopCodon=='TAG' or stopCodon=='TGA' ):
                    try:
                         translationResult2=mRNAsequence.translate(table=11,cds=True)
                    except:
                         translationResult2="Error"
                    cadena="If is not a CDS: "+translationResult1+" . If is a CDS: " + translationResult2
               else:
                    translationResult2="This is not a valid Stop Codon"
                    cadena="If is not a CDS: "+translationResult1+" . If is a CDS: " + translationResult2
               
               results.append(cadena)#2 Translation Result
          else:
               translationResult=DNAsequence.translate(table=codonTable)
               results.append(translationResult)#2 Translation Result
          return(results)


     elif isProt(sequence)==True:
          results.append("This is a nucleotide sequence")#1 type of sequence-Proteins
          results.append(Seq(sequence))#2- this not requires to compute the translation operation
          return(results)
     else:
          results.append("This is not a recognized biological sequence.")#1 type of sequence-NONE
          results.append("This is not a recognized biological sequence.")#2 results of transcription
          return(results)

################################################################################################################
#def sequenceAligner(sequence1,sequence2,mode,subsMat,mscore,mismscore,tiogs,tiegs,tlegs,tlogs,tregs,trogs,qiegs,qiogs,qlegs,qlogs,qregs,qrogs,gapScore):
def sequenceAligner(sequence1,sequence2,mode,subsMat,mscore,mismscore,gapScore):
     sequence1=sequence1.replace('*','')
     sequence2=sequence2.replace('*','')

     if mscore=='':
          mscore=1
     if mismscore=='':
          mismscore=0
          
     if gapScore=='':
          gapScore=0
     

     mscore=float(mscore)

     gapScore=float(gapScore)
     if mscore==0:
          mscore=0.000001
     mismscore=float(mismscore)
     
     results=[]
     sequence1=sequence1.upper()
     sequence2=sequence2.upper()
     sequence1=Seq(sequence1)
     results.append(sequence1)#0-sequence 1
     sequence2=Seq(sequence2)
     results.append(sequence2)#1-sequence 2
     aligner=Align.PairwiseAligner()

     if subsMat!="":
          matrix=Align.substitution_matrices.load(subsMat)
          aligner.substitution_matrix=matrix
     
     aligner.match_score=mscore
     aligner.mismatch_score=mismscore
     
     aligner.gap_score=gapScore

     aligner.mode=mode#global or local
     alignments=aligner.align(sequence1,sequence2)
     results.append(alignments)#2-alignments object
     score=aligner.score(sequence1,sequence2)
     lista=[]
     for alignment in alignments:
          cadena=format_alignment(*alignment,score=score,begin=0,end=alignment.shape[1])
          lista.append(alignment)
     results.append(cadena)#3-line by line alignment
     results.append(lista)#4-lista append con \n
     #print(alignments)
     score=aligner.score(sequence1,sequence2)
     results.append(score)#5-score of alignments  
     results.append(aligner)#6-parammeters of alignments
     results.append(mscore)#7-match score
     results.append(subsMat)#8-substitution_matrix
     results.append(mode)#9-substitution_matrix
     results.append(mismscore)#10-mismatchScore

     results.append(gapScore)#23/11

     return(results)        


def ncbiBlastToolFunction(sequence,type,passD):
     sequence=sequence.replace('*','')
     #print("Running NCBI Blast Tool")
     results=[]
     results.append(sequence)#0 - sequence or number of GI identifier

     results.append(type)#1 - type of search blastp, blastn, ...
     resultHandle=NCBIWWW.qblast(type,passD,sequence)
     blastRecords = NCBIXML.parse(resultHandle)
     blastData=[]
     for record in blastRecords:
          for alignment in record.alignments:
               for hsp in alignment.hsps:
                    blastData.append({
                         'Query ID': record.query_id,
                         'Alignment Title': alignment.title,
                         'Alignment Length': alignment.length,
                         'E-value': hsp.expect,
                         'Score': hsp.score,
                         'Query Sequence': hsp.query,
                         'Alignment Sequence': hsp.sbjct
                         })
     df=pd.DataFrame(blastData)
     df=df.to_html().replace('<td>', '<td align="center">')
     df=df.replace('<th>', '<th align="center">')
     results.append(df)#2 - df
     return(results)


fileBytes=bytes()
def philogeneticTreeFunction(f):
     calculator = DistanceCalculator('identity')
     results=[]
     
     alineamiento=[]

     with open('main/media/upload/' + f.name ,'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)
            fileBytes=chunk
     #print("Archivo subido")

     #print(f.name)

     with open("main/media/upload/" + f.name, "r") as aln:
          alignment = AlignIO.read(aln, "clustal")
          alineamiento.append(alignment)

     results.append(alineamiento)#0 -alineamiento
     dm = calculator.get_distance(alignment)
     results.append(dm)#1- calculadora de alineamientos

     #initialize a DistanceTreeConstructor object based on our distance calculator
     object 
     constructor = DistanceTreeConstructor(calculator)
     #build the tree
     tree = constructor.build_tree(alignment)
     arbolGenerado=Phylo.draw(tree)

     results.append(Phylo.draw(tree))#2
     results.append(tree)#3

     #print(tree)
     buf = BytesIO()
     grafica=Phylo.draw(tree, do_show=False)
     plt.savefig(buf,format='png')
     image_base64 = base64.b64encode(buf.getvalue()).decode('utf-8').replace('\n', '')
     buf.close()
     results.append(image_base64)#4

     os.remove("main/media/upload/"+f.name)
     return(results)


def pdbAnalysis(f):
     results=[]
     model_list=[]
     chain_list=[]
     residues_list=[]
     atom_name=[]
     atom_vectors=[]
     unified_list=[]

     pdb_parser = PDB.PDBParser(QUIET=True)
     
     
     with open('main/media/upload/' + f.name ,'wb+') as destination:
          for chunk in f.chunks():
               destination.write(chunk)
               fileBytes=chunk
               
     structure = pdb_parser.get_structure("PDB_File", 'main/media/upload/'+f.name)
     ########################################## BIOPANDAS ANALYSIS ##########################################
     structure_biopandas=PandasPdb().read_pdb('main/media/upload/'+f.name)
     atom_df=structure_biopandas.df['ATOM']

     het_df = structure_biopandas.df['HETATM']

     buf = BytesIO()

     atom_df['b_factor'].plot(kind='hist')     
     plt.xlabel="B-Factor"
     plt.title="B-Factor"
          
     plt.savefig(buf, format='jpg')
     betaFactorPlot = base64.b64encode(buf.getvalue()).decode('utf-8').replace('\n', '')
     buf.close()
     plt.close()

     buf1 = BytesIO()
     plt.clf()
     atom_df['element_symbol'].value_counts().plot(kind='bar')
     
     plt.ylabel('Count')
     
     plt.savefig(buf1,format='jpg')
     elementDistributionPlot=base64.b64encode(buf1.getvalue()).decode('utf-8').replace('\n', '')
     buf1.close()
     plt.close()
     
     buf2 = BytesIO()
     plt.clf()
     
     atom_df['atom_name'].value_counts().plot(kind='bar')
     plt.ylabel('Count')
     plt.savefig(buf2,format='jpg')
     
     atomDistributionPlot=base64.b64encode(buf2.getvalue()).decode('utf-8').replace('\n', '')
     
     buf2.close()
     plt.close()
     ########################################################################################################

     structure_name=structure.header["name"]
     results.append(structure_name)#0-Nombre de la estructura
     depositation_date=structure.header['deposition_date']
     results.append(depositation_date)#1-depositation_date
     structure_release_date=structure.header["release_date"]
     results.append(structure_release_date)#2-Fecha de salida de la estructura
     structure_resolution = structure.header["resolution"]
     results.append(structure_resolution)#3-Resolución de la estructura
     structure_keywords = structure.header["keywords"]
     results.append(structure_keywords)#4-keywords de la estructura

     structure_method=structure.header["structure_method"]
     results.append(structure_method)#5-structure_method
     structure_reference=structure.header['structure_reference']
     results.append(structure_reference)#6-structure_reference
     journal_reference=structure.header['journal_reference']
     results.append(journal_reference)#7-journal_reference
     author=structure.header['author']
     results.append(author)#8-author
     compound=structure.header['compound']
     results.append(compound)#9-compound
     source=structure.header['source']
     results.append(source)#10-source
     has_missing_residues=structure.header["has_missing_residues"]
     results.append(has_missing_residues)#11-has_missing_resiudes
     missing_residues=structure.header["missing_residues"]
     results.append(missing_residues)#12-missing_residues
     journal=structure.header['journal']
     results.append(journal)#13-journal

     for model in structure:
          for chain in model:
               for residue in chain:
                    if 'NGT' in residue.id:
                         cadenaGlycosilation="Glycosylation Found: Model"+ model.id +  "Chain" + chain.id +", Residue "+ residue.id
                    else:
                         cadenaGlycosilation="No information found"
     
     results.append(cadenaGlycosilation)#14-Glycosilation

     
     for model in structure.get_models():
          model_list.append(model)
          for chain in model.get_chains():
               chain_list.append(chain)
               for residue in chain.get_residues():
                    residues_list.append(residue)
                    for atom in residue.get_atoms():
                         atom_name.append(atom)
                         atom_vectors.append(atom.get_vector())
                    break
               break
          break

     for i in range(len(atom_name)):
          string_uni=atom_name[i] , atom_vectors[i]
          unified_list.append(string_uni)

     
     results.append(model_list)#15-model_list
     results.append(chain_list)#16-chains
     results.append(residues_list)#17-residues_list
     results.append(unified_list)#18-atom_name & atom_vectors

     results.append(atom_df)#19-dataframe of atoms
     results.append(het_df)#20- hetero-atom
     results.append(betaFactorPlot)#21 beta-factor plot
     results.append(elementDistributionPlot)#22 elementDistributionPlot
     results.append(atomDistributionPlot)#23 atomDistributionPlot



     os.remove("main/media/upload/"+f.name)
     return results

def loadDataProtFungi():
     lista=[]
     connect_string = "mongodb://localhost:27017/"
     my_client = pymongo.MongoClient(connect_string)
     mongo_db=my_client.fungiRegExAlgal

     collection=mongo_db.protFungi

     for document in collection.find():
          lista.append(document)
     


     return(lista)
     
def loadDataProtAlgal():
     lista=[]
     connect_string = "mongodb://localhost:27017/"
     my_client = pymongo.MongoClient(connect_string)
     mongo_db=my_client.fungiRegExAlgal

     collection=mongo_db.protAlgal

     for document in collection.find():
          lista.append(document)

     
     return(lista)


def loadDataGeneFungi():
     lista=[]
     connect_string = "mongodb://localhost:27017/"
     my_client = pymongo.MongoClient(connect_string)
     mongo_db=my_client.fungiRegExAlgal

     collection=mongo_db.geneFungi

     for document in collection.find():
          lista.append(document)

     
     return(lista)
     
def loadDataGeneAlgal():
     lista=[]
     connect_string = "mongodb://localhost:27017/"
     my_client = pymongo.MongoClient(connect_string)
     mongo_db=my_client.fungiRegExAlgal

     collection=mongo_db.geneAlgal

     for document in collection.find():
          lista.append(document)

     
     return(lista)

def motifAnalysisFunction(sequences):
     list=[]
     results=[]
     try:
          results.append(sequences)#0 sequences list
          sequences=sequences.replace(" ","")
          sequences=sequences.split(",")
          for element in sequences:
               element=element.replace(" ","")
               list.append(Seq(element))
               
          m=motifs.create(list)
               
          conteos=m.counts
          results.append(conteos)#1-matriz de conteos
          #print(conteos)
          consenso=m.consensus
          results.append(consenso)#2-consenso
          degenerado_consenso=m.degenerate_consensus
          results.append(degenerado_consenso)#3-consenso degenerado
          reversa_complementaria=m.reverse_complement()
          results.append(reversa_complementaria)#4-reversa complementaria
          consensus_reverse_complement=reversa_complementaria.consensus
          results.append(consensus_reverse_complement)#5-consenso de reversa complementaria
          return(results)
     except:
          i=0
          for i in range(5):
               results.append("Error, check your motifs and documentation.")

          return(results)


def regularExpressionModule(database,regularExpression):
     matches=[]
     numMatches=[]
     maxi=0
     if database=="fungiProteome":
          
          data=loadDataProtFungi()
          for element in data:
               search=re.findall(regularExpression,element["seq"])
               if(search):
                    tam=len(search)
                    if tam>maxi:
                         maxi=tam
                    element["numMatches"]=tam
                    matches.append(element)
                    numMatches.append(tam)

     if database=="fungiGenome":
          
          data=loadDataGeneFungi()
          for element in data:
               search=re.findall(regularExpression,element["seq"])
               if(search):
                    tam=len(search)
                    if tam>maxi:
                         maxi=tam
                    element["numMatches"]=tam
                    matches.append(element)
                    numMatches.append(tam)
          
     if database=="algalGenome":
          
          data=loadDataGeneAlgal()
          for element in data:
               search=re.findall(regularExpression,element["seq"])
               if(search):
                    tam=len(search)
                    if tam>maxi:
                         maxi=tam
                    element["numMatches"]=tam
                    matches.append(element)
                    numMatches.append(tam)
     if database=="algalProteome":
          
          data=loadDataProtAlgal()
          for element in data:
               search=re.findall(regularExpression,element["seq"])
               if(search):
                    tam=len(search)
                    if tam>maxi:
                         maxi=tam
                    element["numMatches"]=tam
                    matches.append(element)
                    numMatches.append(tam)
     
     return(matches,maxi)