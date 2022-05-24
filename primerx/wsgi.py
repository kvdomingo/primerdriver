from .app import app as application
from .log import logger


@logger.catch
def main():
    application.run()


if __name__ == "__main__":
    main()
