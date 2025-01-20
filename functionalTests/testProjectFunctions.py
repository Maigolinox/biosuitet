from selenium import webdriver
from django.contrib.staticfiles.testing import StaticLiveServerTestCase
from selenium.webdriver.chrome.service import Service
from django.urls import reverse
import time


class testProjectFunctions(StaticLiveServerTestCase):

    def setUp(self):
        self.browser =webdriver.Chrome('functionalTests/chromedriver.exe')

    def tearDown(self):
        self.browser.close()


    def testFunction1(self):
        self.browser.get(self.live_server_url)
        time.sleep(20)