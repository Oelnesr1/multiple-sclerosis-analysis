import selenium
from selenium import webdriver
from selenium import *
from selenium.webdriver.opera.options import Options
import re
import time

DRIVER_PATH = '/Users/omarelnesr/operadriver_mac64/operadriver'
url = 'https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA727413&o=acc_s%3Aa'
options = Options()
options.headless = True
options.add_argument("--window-size=1920,1200")

driver = webdriver.Opera(options=options, executable_path=DRIVER_PATH)


SRRs = []
for i in range(1):
    SRRs.append('SRR' + str(14409206 + i))
# driver.get(url)
# time.sleep(10)
# html = driver.page_source
# tab = '\t'
print(SRRs)
# for SRR in SRRs:
#     print(SRR)

