from json2html import *
import re
import json
from django.shortcuts import render, redirect

from django.utils.decorators import method_decorator
from django.contrib.auth.decorators import login_required
from allauth.socialaccount.models import SocialAccount

from django.http import HttpResponseRedirect#redireccionar a la misma pagina

from django.urls import reverse#redireccionar a dashboard

import os

from .functions import biotraToolFunction, proteinAnalysisF, transcriptionFunction, backTranscriptionTool, translationTool, sequenceAligner, ncbiBlastToolFunction, philogeneticTreeFunction, pdbAnalysis, loadDataProtFungi, loadDataGeneFungi, loadDataGeneAlgal, loadDataProtAlgal, motifAnalysisFunction, regularExpressionModule

from .forms import UploadFileForm




# Create your views here.
def home(request):
    if request.method=="POST":
        pass
    else:
        return render(request,'home.html')
    
def dashboard(request):
    if request.user.is_authenticated:
        return render(request,'dashboard.html')
    else:
        return render(request,'forbiden.html')

def biotratool(request):
    sendGeneralData=[]
    if request.user.is_authenticated:
        if request.method=="POST":
            sequence=request.POST.get('sequence')
            result=biotraToolFunction(sequence)
            sendGeneralData.append(result[3]*100)#A
            sendGeneralData.append(result[4]*100)#G
            sendGeneralData.append(result[5]*100)#C
            sendGeneralData.append(result[6]*100)#T
            return render(request,'propertiesTool.html',context={'result':result,'sequence':sequence,'sendGeneralData':sendGeneralData})
        else:
            return render(request,'propertiesTool.html')
    else:
        return render(request,'forbiden.html')
    
def proteinAnalysisView(request):
    sendGeneralData=[]
    if request.user.is_authenticated:
        if request.method=="POST":
            sequence=request.POST.get('sequence')
            pHValue=request.POST.get('pHValue')
            try:
                resultsProt=proteinAnalysisF(sequence,pHValue)
            except:
                resultsProt="Error"
            
            def isDNA(seq, alphabet='dna'):
             alphabets = {'dna': re.compile('^[acgtn]*$', re.I), 'protein': re.compile('^[acdefghiklmnpqrstvwy]*$', re.I)}
             if alphabets[alphabet].search(seq) is not None:
                  return True
             else:
                return False
            
            if isDNA(sequence)==True:
                sendGeneralData=0
                keys=0
                values=0
            else:
                try:
                    sendGeneralData={x:y for x,y in resultsProt[1].items() if y!=0}
                    keys=list(sendGeneralData.keys())
                    values=list(sendGeneralData.values())
                except:
                    sendGeneralData="Error"
                    keys="Error"
                    values="Error"
                
            #sendGeneralData=list(resultsProt[1].values())
            
            return render(request,'propertiesToolProt.html',context={'resultsProt':resultsProt,'sequence':sequence,'sendGeneralData':sendGeneralData,'pHValue':pHValue,'keys':keys,'values':values})
        else:
            return render(request,'propertiesToolProt.html')
    else:
        return render(request,'forbiden.html')
    

def transcriptionToolView(request):
    if request.user.is_authenticated:
        if request.method=="POST":
            DNAsequence=request.POST.get('DNAsequence')
            
            resultsTranscriptionFunction=transcriptionFunction(DNAsequence)
            return render(request,'transcriptionTool.html',context={'resultsTranscriptionFunction':resultsTranscriptionFunction})
        else:
            return render(request,'transcriptionTool.html')
    else:
        return render(request,'forbiden.html')
    
def backTranscriptionToolView(request):
    if request.user.is_authenticated:
        if request.method=="POST":
            mRNAsequence=request.POST.get('mRNASequence')
            
            resultsBackTranscriptionFunction=backTranscriptionTool(mRNAsequence)
            return render(request,'backTranscriptionTool.html',context={'resultsBackTranscriptionFunction':resultsBackTranscriptionFunction})
        else:
            return render(request,'backTranscriptionTool.html')
    else:
        return render(request,'forbiden.html')


def translationToolView(request):
    if request.user.is_authenticated:
        if request.method=="POST":
            sequence=request.POST.get('sequence')
            codonTable=request.POST.get('codonTable')
            
            codonTable=int(codonTable)
            resultsTranslationTool=translationTool(sequence,codonTable)
            return render(request,"translationTool.html",context={'resultsTranslationTool':resultsTranslationTool,'codonTable':codonTable})
        else:
            return render(request,"translationTool.html")
        
    else:
        return render(request,'forbiden.html')
    
def sequenceAlignerToolView(request):
    if request.user.is_authenticated:
        if request.method=="POST":
            sequence1=request.POST.get('sequence1')
            sequence2=request.POST.get('sequence2')
            subsMat=request.POST.get('subsMat')
            mode=request.POST.get('mode')
            mscore=request.POST.get('mscore')
            mismscore=request.POST.get('mismscore')
            
            gapScore=request.POST.get('gapScore')
            

            resultsSequenceAligner=sequenceAligner(sequence1,sequence2,mode,subsMat,mscore,mismscore,gapScore)
            return render(request,"sequenceAlignerTool.html",context={'resultsSequenceAligner':resultsSequenceAligner})
        else:
            return render(request,"sequenceAlignerTool.html")
        
    else:
        return render(request,'forbiden.html')
    


def blastToolView(request):
    if request.user.is_authenticated:
        if request.method=="POST":
            sequence=request.POST.get('sequenceGI')
            #print(sequence)
            passD="nt"
            
            type=request.POST.get('typeB')
            
            resultsNCBI=ncbiBlastToolFunction(sequence,type,passD)

            return render(request,"ncbiBlastTool.html",context={'resultsNCBI':resultsNCBI})
        else:
            return render(request,"ncbiBlastTool.html")
    else:
        return render(request,'forbiden.html')
    
    
def phyloGeneticTreesToolView(request):
    if request.user.is_authenticated:
        
        if request.method=="POST":
            form=UploadFileForm(request.POST, request.FILES)
            if form.is_valid():
                #print("Formulario válido")
                #print(form)
                try:
                    results=philogeneticTreeFunction(request.FILES['file'])
                except:
                    results=["Error creating philogenetic tree","Error"]
                #print(results)
            else:
                #print("Formulario Invalido")
                
                results=["Error creating philogenetic tree","Error"]
                
            
            return render(request,"phyloGeneticTreesTool.html", context={'form':form,'results':results})
        else:
            form=UploadFileForm()
            return render(request,"phyloGeneticTreesTool.html",context={'form':form})
    else:
        return render(request,'forbiden.html')
    

def pdbViewerToolView(request):
    if request.user.is_authenticated:
        if request.method=="POST":
            return render(request,"pdbViewerTool.html")
        else:
            return render(request,"pdbViewerTool.html")
    else:
        return render(request,'forbiden.html')

def pdbAnalysisToolView(request):
    if request.user.is_authenticated:
        if request.method=="POST":
            form=UploadFileForm(request.POST, request.FILES)
            if form.is_valid():
                #print("Formulario válido")
                #results=pdbAnalysis(request.FILES['file'])

                try:
                    results=pdbAnalysis(request.FILES['file'])
                except:
                    #results=pdbAnalysis(request.FILES['file'])
                    results=["Error analyzing PDB File","Error"]
            else:
                #print("Formulario Invalido")
                
                results=["Error analyzing PDB file","Error","Error"]
            return render(request,"pdbAnalysisTool.html",context={'form':form,'results':results})
        else:
            form=UploadFileForm()
            return render(request,"pdbAnalysisTool.html",context={'form':form})
    else:
        return render(request,'forbiden.html')


def fungiProtView(request):
    lista=loadDataProtFungi()
    specieList=[]
    IDList=[]
    seqList=[]
    i=0
    
    
    # for i in range(len(lista)):
    #     y=json.loads(lista[i])
    #     specieList.append(y["specie"])
    #     IDList.append(y["ID"])
    #     seqList.append(y["seq"])

    # lista=[]
    #print(type(lista))
    if request.user.is_authenticated:
        if request.method=="POST":
            
            
            return render(request,"fungiProtDB.html")
        else:
            
            
            return render(request,"fungiProtDB.html",context={'specieList':specieList,'IDList':IDList,'seqList':seqList,'lista':lista})
    else:
        return render(request,'forbiden.html')
    
def fungiGeneView(request):
    lista=loadDataGeneFungi()
    if request.user.is_authenticated:
        
        
        return render(request,"fungiGeneTool.html",context={'lista':lista})
    else:
        return render(request,'forbiden.html')
    
def algalProtView(request):
    lista=loadDataProtAlgal()
    
    if request.user.is_authenticated:
        if request.method=="POST":
            
            

            return render(request,"algalProtDB.html",context={'lista':lista})
        else:
            return render(request,"algalProtDB.html",context={'lista':lista})
    else:
        return render(request,'forbiden.html')
    
def algalGeneView(request):
    lista=loadDataGeneAlgal()
    if request.user.is_authenticated:
        if request.method=="POST":
            return render(request,"algalGeneTool.html",context={'lista':lista})
        else:
            return render(request,"algalGeneTool.html",context={'lista':lista})
    else:
        return render(request,'forbiden.html')



def motifsToolView(request):
    if request.user.is_authenticated:
        if request.method=="POST":
            
            motifSequences=request.POST.get('motifSequences')
            #print(motifSequences)
            results=motifAnalysisFunction(motifSequences)
            return render(request,"motifsTool.html",context={'results':results})
        else:
            return render(request,"motifsTool.html")
    else:
        return render(request,'forbiden.html')
    


def fungiRegExView(request):
    if request.user.is_authenticated:
        if request.method=="POST":
            database=request.POST.get('database')
            regularExpression=request.POST.get('regularExpression')
            results,maxMatches=regularExpressionModule(database,regularExpression)
            return render(request,"fungiRegExTool.html",context={'results':results,'database':database,'regularExpression':regularExpression,'maxMatches':maxMatches})
        else:
            return render(request,"fungiRegExTool.html")
    else:
        return render(request,'forbiden.html')
    

def privacyPolicyView(request):
    return render(request,"privacyPolicy.html")