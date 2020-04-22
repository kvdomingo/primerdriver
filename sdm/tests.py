from selenium import webdriver
from percy import percySnapshot
from django.test import TestCase
from django.conf import settings


class LandingTestCase(TestCase):
    browser = webdriver.Firefox()
    browser.get('http://localhost:8000')
    browser.implicitly_wait(10)
    browser.execute_script('window.scrollTo(100, Y)')
    percySnapshot(browser=browser, name='homepage')
    browser.close()
