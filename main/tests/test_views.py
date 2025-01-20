from django.test import TestCase, Client
from django.urls import reverse
import json


class TestViews(TestCase):

    def setUp(self):
        self.client=Client()
        self.homeURL=reverse('homepage')
        self.dashboardURL=reverse('dashboard')
        self.biotratoolURL=reverse('biotratool')
        self.bioprotoolURL=reverse('bioprotool')
        self.transcriptionToolURL=reverse('transcriptionTool')
        self.backTranscriptionToolURL=reverse('backTranscriptionTool')
        self.translationToolURL=reverse('translationTool')
        self.sequenceAlignerToolURL=reverse('sequenceAlignerTool')
        self.ncbiBlastToolURL=reverse('ncbiBlastTool')
        self.fungiProtDBURL=reverse('fungiProtDB')
        self.fungiGenoDBURL=reverse('fungiGenoDB')
        self.algalProtDBURL=reverse('algalProtDB')
        self.algalGenoDBURL=reverse('algalGenoDB')
        self.phyloTreesToolURL=reverse('phyloTreesTool')
        self.pdbViewerToolURL=reverse('pdbViewerTool')
        self.pdbAnalysisToolURL=reverse('pdbAnalysisTool')
        self.regExToolURL=reverse('regExTool')
        self.motifsToolURL=reverse('motifsTool')
        self.privacyPolicyURL=reverse('privacyPolicy')




    def testHome_GET(self):
        client = Client()
        response = client.get(self.homeURL)
        self.assertEquals(response.status_code,200)        
        self.assertTemplateUsed(response,'home.html')

    def testDashboard_GET(self):
        client = Client()
        response = client.get(self.dashboardURL)
        self.assertEquals(response.status_code,200)        
        self.assertTemplateUsed(response,'forbiden.html')

    def biotratool_GET(self):
        client = Client()
        response = client.get(self.biotratoolURL)
        self.assertEquals(response.status_code,200)        
        self.assertTemplateUsed(response,'propertiesTool.html')

    def proteinAnalysisView_GET(self):
        client = Client()
        response = client.get(self.bioprotoolURL)
        self.assertEquals(response.status_code,200)        
        self.assertTemplateUsed(response,'propertiesToolProt.html')

    def backTranscriptionToolView_GET(self):
        client = Client()
        response = client.get(self.backTranscriptionToolURL)
        self.assertEquals(response.status_code,200)        
        self.assertTemplateUsed(response,'backTranscriptionTool.html')

    def translationToolView_GET(self):
        client = Client()
        response = client.get(self.translationToolURL)
        self.assertEquals(response.status_code,200)        
        self.assertTemplateUsed(response,'translationTool.html')


    def sequenceAlignerToolView_GET(self):
        client = Client()
        response = client.get(self.translationToolURL)
        self.assertEquals(response.status_code,200)        
        self.assertTemplateUsed(response,'sequenceAlignerTool.html')
    
    def blastToolView_GET(self):
        client = Client()
        response = client.get(self.ncbiBlastToolURL)
        self.assertEquals(response.status_code,200)        
        self.assertTemplateUsed(response,'ncbiBlastTool.html')

    def fungiProtView_GET(self):
        client = Client()
        response = client.get(self.fungiProtDBURL)
        self.assertEquals(response.status_code,200)        
        self.assertTemplateUsed(response,'fungiProtDB.html')

    def fungiGeneView_GET(self):
        client = Client()
        response = client.get(self.fungiGenoDBURL)
        self.assertEquals(response.status_code,200)        
        self.assertTemplateUsed(response,'fungiGeneDB.html')

    def algalProtView_GET(self):
        client = Client()
        response = client.get(self.algalProtDBURL)
        self.assertEquals(response.status_code,200)        
        self.assertTemplateUsed(response,'algalProtDB.html')

    def algalGeneView_GET(self):
        client = Client()
        response = client.get(self.algalGenoDBURL)
        self.assertEquals(response.status_code,200)        
        self.assertTemplateUsed(response,'algalGeneTool.html')

    def phyloGeneticTreesToolView_GET(self):
        client = Client()
        response = client.get(self.phyloTreesToolURL)
        self.assertEquals(response.status_code,200)        
        self.assertTemplateUsed(response,'phyloGeneticTreesTool.html')

    def pdbViewerToolView_GET(self):
        client = Client()
        response = client.get(self.pdbViewerToolURL)
        self.assertEquals(response.status_code,200)        
        self.assertTemplateUsed(response,'pdbViewerTool.html')

    def pdbAnalysisToolView_GET(self):
        client = Client()
        response = client.get(self.pdbAnalysisToolURL)
        self.assertEquals(response.status_code,200)        
        self.assertTemplateUsed(response,'pdbAnalysisTool.html')

    def fungiRegExView_GET(self):
        client = Client()
        response = client.get(self.phyloTreesToolURL)
        self.assertEquals(response.status_code,200)        
        self.assertTemplateUsed(response,'fungiRegExTool.html')

    def motifsToolView_GET(self):
        client = Client()
        response = client.get(self.phyloTreesToolURL)
        self.assertEquals(response.status_code,200)        
        self.assertTemplateUsed(response,'motifsTool.html')

    def privacyPolicyView(self):
        client=Client()
        response=client.get(self.privacyPolicyURL)
        self.assertEquals(response.status_code,200)
        self.assertTemplateUsed(response,'privacyPolicy.html')