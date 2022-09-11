wsgi_app = "primerx.wsgi"

worker_class = "gthread"
workers = 1
threads = 2

errorlog = "-"
accesslog = "-"
capture_output = True
loglevel = "error"
