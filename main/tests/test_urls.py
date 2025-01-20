from django.test import SimpleTestCase
from django.urls import reverse, resolve
from main.views import home,dashboard, biotratool,proteinAnalysisView, transcriptionToolView,backTranscriptionToolView,translationToolView, sequenceAlignerToolView, blastToolView,fungiProtView,fungiGeneView,algalProtView, algalGeneView,phyloGeneticTreesToolView,pdbViewerToolView,pdbAnalysisToolView,fungiRegExView,motifsToolView,privacyPolicyView

class testUrls(SimpleTestCase):
    def testHomepageUrlIsResolved(self):
        url=reverse('homepage')
        print(resolve(url))
        self.assertEquals(resolve(url).func,home)

    def testDashboardUrlIsResolved(self):
        url=reverse('dashboard')
        print(resolve(url))
        self.assertEquals(resolve(url).func,dashboard)

    def testBioTraToolUrlIsResolved(self):
        url=reverse('biotratool')
        print(resolve(url))
        self.assertEquals(resolve(url).func,biotratool)

    def testBioproToolUrlIsResolved(self):
        url=reverse('bioprotool')
        print(resolve(url))
        self.assertEquals(resolve(url).func,proteinAnalysisView)

    def testTranscriptionUrlIsResolved(self):
        url=reverse('transcriptionTool')
        print(resolve(url))
        self.assertEquals(resolve(url).func,transcriptionToolView)

    def testBackTranscriptionUrlIsResolved(self):
        url=reverse('backTranscriptionTool')
        print(resolve(url))
        self.assertEquals(resolve(url).func,backTranscriptionToolView)

    def testTranslationToolUrlIsResolved(self):
        url=reverse('translationTool')
        print(resolve(url))
        self.assertEquals(resolve(url).func,translationToolView)

    def testSequenceAlignerToolUrlIsResolved(self):
        url=reverse('sequenceAlignerTool')
        print(resolve(url))
        self.assertEquals(resolve(url).func,sequenceAlignerToolView)

    def testBlastToolUrlIsResolved(self):
        url=reverse('ncbiBlastTool')
        print(resolve(url))
        self.assertEquals(resolve(url).func,blastToolView)

    def testFungiProtToolUrlIsResolved(self):
        url=reverse('fungiProtDB')
        print(resolve(url))
        self.assertEquals(resolve(url).func,fungiProtView)

    def testFungiGeneToolUrlIsResolved(self):
        url=reverse('fungiGenoDB')
        print(resolve(url))
        self.assertEquals(resolve(url).func,fungiGeneView)

    def testAlgalProtUrlIsResolved(self):
        url=reverse('algalProtDB')
        print(resolve(url))
        self.assertEquals(resolve(url).func,algalProtView)

    def testAlgalGeneUrlIsResolved(self):
        url=reverse('algalGenoDB')
        print(resolve(url))
        self.assertEquals(resolve(url).func,algalGeneView)

    def testPhyloUrlIsResolved(self):
        url=reverse('phyloTreesTool')
        print(resolve(url))
        self.assertEquals(resolve(url).func,phyloGeneticTreesToolView)

    def testPdbViewerUrlIsResolved(self):
        url=reverse('pdbViewerTool')
        print(resolve(url))
        self.assertEquals(resolve(url).func,pdbViewerToolView)

    def testPdbAnalysisToolUrlIsResolved(self):
        url=reverse('pdbAnalysisTool')
        print(resolve(url))
        self.assertEquals(resolve(url).func,pdbAnalysisToolView)

    def testRegExToolUrlIsResolved(self):
        url=reverse('regExTool')
        print(resolve(url))
        self.assertEquals(resolve(url).func,fungiRegExView)

    def testMotifsToolUrlIsResolved(self):
        url=reverse('motifsTool')
        print(resolve(url))
        self.assertEquals(resolve(url).func,motifsToolView)

    def testPrivacyPolicyUrlIsResolved(self):
        url=reverse('privacyPolicy')
        print(resolve(url))
        self.assertEquals(resolve(url).func,privacyPolicyView)