from selenium import webdriver
from selenium.webdriver.firefox.options import Options
from percy import percySnapshot
from django.test import TestCase
from django.conf import settings


class LandingTestCase(TestCase):
    options = Options()
    options.headless = True
    browser = webdriver.Firefox(options=options)
    browser.get('http://0.0.0.0:8000')
    browser.implicitly_wait(10)
    browser.execute_script('window.scrollTo(100, Y)')
    percySnapshot(browser=browser, name='homepage')
    browser.close()
