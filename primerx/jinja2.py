import os
from django.templatetags.static import static
from django.urls import reverse
from django.conf import settings
from webpack_loader.templatetags.webpack_loader import render_bundle
from pdcli.version import __version__
from jinja2 import Environment
from datetime import datetime
from dotenv import load_dotenv


load_dotenv()

def environment(**options):
    env = Environment(**options)
    env.globals.update({
        'now': datetime.now(),
        'program_version': __version__,
        'render_bundle': render_bundle,
        'settings': settings,
        'static': static,
        'url': reverse,
        'web_version': f'(web {os.environ["HEROKU_RELEASE_VERSION"]})' if settings.ON_HEROKU else '',
    })
    return env
