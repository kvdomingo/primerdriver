from django.conf import settings
from django.test import TestCase, LiveServerTestCase
from django.contrib.staticfiles.testing import StaticLiveServerTestCase
from selenium.webdriver.firefox.webdriver import WebDriver
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.support.ui import Select
from percy import percySnapshot


PERCY_BS_WIDTHS = [576, 768, 992, 1200]

class WebTestCase(StaticLiveServerTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        if settings.ON_CI:
            options = Options()
            options.headless = True
            cls.selenium = WebDriver(options=options)
        else:
            cls.selenium = WebDriver()
        cls.selenium.implicitly_wait(10)

    @classmethod
    def tearDownClass(cls):
        cls.selenium.quit()
        super().tearDownClass()

    def test_landing_view(self):
        self.selenium.get(self.live_server_url)
        self.selenium.execute_script('window.scrollTo(0, 100)')
        percySnapshot(browser=self.selenium, name='homepage', widths=PERCY_BS_WIDTHS)

    def test_characterize_view(self):
        self.selenium.get(self.live_server_url)
        goto_char = self.selenium.find_element_by_id('characterize')
        goto_char.click()
        percySnapshot(browser=self.selenium, name='characterize', widths=PERCY_BS_WIDTHS)

    def test_dna_view(self):
        self.selenium.get(self.live_server_url)
        goto_char = self.selenium.find_element_by_id('dna-based')
        goto_char.click()
        percySnapshot(browser=self.selenium, name='dna', widths=PERCY_BS_WIDTHS)

    def test_protein_view(self):
        self.selenium.get(self.live_server_url)
        goto_char = self.selenium.find_element_by_id('protein-based')
        goto_char.click()
        percySnapshot(browser=self.selenium, name='protein', widths=PERCY_BS_WIDTHS)

    def test_characterize_submit_disabled_on_characterize_landing(self):
        self.selenium.get(self.live_server_url)
        goto_char = self.selenium.find_element_by_id('characterize')
        goto_char.click()
        submit_button = self.selenium.find_element_by_id('submit')
        submit_disabled = submit_button.get_attribute('disabled')
        self.assertEqual(submit_disabled, 'true')

    def test_characterize_submit_disabled_on_type_sequence(self):
        self.selenium.get(self.live_server_url)
        goto_char = self.selenium.find_element_by_id('characterize')
        goto_char.click()
        textarea = self.selenium.find_element_by_id('sequence')
        textarea.send_keys('GATTACA')
        submit_button = self.selenium.find_element_by_id('submit')
        submit_disabled = submit_button.get_attribute('disabled')
        self.assertEqual(submit_disabled, 'true')

    def test_characterize_submit_disabled_on_set_mismatch(self):
        self.selenium.get(self.live_server_url)
        goto_char = self.selenium.find_element_by_id('characterize')
        goto_char.click()
        textarea = self.selenium.find_element_by_id('sequence')
        textarea.send_keys('GATTACA')
        mismatch = self.selenium.find_element_by_id('mismatched_bases')
        mismatch.send_keys('1')
        submit_button = self.selenium.find_element_by_id('submit')
        submit_disabled = submit_button.get_attribute('disabled')
        self.assertEqual(submit_disabled, 'true')

    def test_characterize_submit_enabled_on_set_mutation(self):
        self.selenium.get(self.live_server_url)
        goto_char = self.selenium.find_element_by_id('characterize')
        goto_char.click()
        textarea = self.selenium.find_element_by_id('sequence')
        textarea.send_keys('GATTACA')
        mismatch = self.selenium.find_element_by_id('mismatched_bases')
        mismatch.send_keys('1')
        mutation = Select(self.selenium.find_element_by_id('mutation_type'))
        mutation.select_by_visible_text('Substitution')
        submit_button = self.selenium.find_element_by_id('submit')
        submit_disabled = submit_button.get_attribute('disabled')
        self.assertEqual(submit_disabled, None)

    def test_characterize_submit_disabled_on_form_reset(self):
        self.selenium.get(self.live_server_url)
        goto_char = self.selenium.find_element_by_id('characterize')
        goto_char.click()
        textarea = self.selenium.find_element_by_id('sequence')
        textarea.send_keys('GATTACA')
        mismatch = self.selenium.find_element_by_id('mismatched_bases')
        mismatch.send_keys('1')
        mutation = Select(self.selenium.find_element_by_id('mutation_type'))
        mutation.select_by_visible_text('Substitution')
        reset_button = self.selenium.find_element_by_id('reset')
        reset_button.click()
        submit_button = self.selenium.find_element_by_id('submit')
        submit_disabled = submit_button.get_attribute('disabled')
        self.assertEqual(submit_disabled, 'true')

    def test_characterize_press_back(self):
        self.selenium.get(self.live_server_url)
        goto_char = self.selenium.find_element_by_id('characterize')
        goto_char.click()
        back_button = self.selenium.find_element_by_id('back')
        back_button.click()
        submit_buttons = self.selenium.find_elements_by_id('submit')
        self.assertEqual(len(submit_buttons), 0)

    def test_characterize_submit(self):
        self.selenium.get(self.live_server_url)
        goto_char = self.selenium.find_element_by_id('characterize')
        goto_char.click()
        textarea = self.selenium.find_element_by_id('sequence')
        textarea.send_keys('CGATCGTACGGACGCAGCTCGTAGCTACGATCGATCGATCGATCGTACGTACGTACGATCGTACGATCGATCGTACG')
        mismatch = self.selenium.find_element_by_id('mismatched_bases')
        mismatch.send_keys('1')
        mutation = Select(self.selenium.find_element_by_id('mutation_type'))
        mutation.select_by_visible_text('Substitution')
        submit = self.selenium.find_element_by_id('submit')
        submit.click()
        table_rows = self.selenium.find_elements_by_tag_name('td')
        self.assertEqual(len(table_rows), 14)
