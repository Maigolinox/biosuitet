from django.urls import path
from . import views
from django.contrib.auth import views as auth_views
from django.contrib.staticfiles.urls import staticfiles_urlpatterns
from django.conf import settings
from django.conf.urls.static import static
from . import views
from . import templates

urlpatterns=[
    path('',views.home,name="homepage"),
    path('dashboard/',views.dashboard,name="dashboard"),
    path('bioprotool/',views.biotratool,name="biotratool"),
    path('bioprotoolProt/',views.proteinAnalysisView,name="bioprotool"),
    path('transcriptionTool/',views.transcriptionToolView,name="transcriptionTool"),
    path('backTranscriptionTool/',views.backTranscriptionToolView,name="backTranscriptionTool"),
    path('translationTool/',views.translationToolView,name="translationTool"),
    path('sequenceAlignerTool/',views.sequenceAlignerToolView,name="sequenceAlignerTool"),
    path('ncbiBlastTool/',views.blastToolView,name="ncbiBlastTool"),
    path('fungiProtDB/',views.fungiProtView,name="fungiProtDB"),
    path('fungiGenoDB/',views.fungiGeneView,name="fungiGenoDB"),
    path('algalProtDB/',views.algalProtView,name="algalProtDB"),
    path('algalGenoDB/',views.algalGeneView,name="algalGenoDB"),
    path('phyloTreesTool/',views.phyloGeneticTreesToolView,name="phyloTreesTool"),
    path('pdbViewerTool/',views.pdbViewerToolView,name="pdbViewerTool"),
    path('pdbAnalysisTool/',views.pdbAnalysisToolView,name="pdbAnalysisTool"),
    path('regExTool/',views.fungiRegExView,name="regExTool"),
    path('motifsTool/',views.motifsToolView,name="motifsTool"),
    path('privacyPolicy/',views.privacyPolicyView,name="privacyPolicy"),




]