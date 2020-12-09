import mechanicalsoup
from bs4 import BeautifulSoup

class SMILEStrans(object):
    def __init__(self, smiles = ''):
        self.br = mechanicalsoup.StatefulBrowser(user_agent='Mozilla/5.0')
        self.br.open('https://cactus.nci.nih.gov/translate/index.html#Form',verify=False)
        self.SMILES = smiles
        self.brtxt = None
        self.trans = ''
        self.resp = None
        self.result = ''

    def setSMILES(self, smiles):
        if isinstance(smiles, str):
            self.__init__()
            self.SMILES = smiles
            print("SMILES: %s" %(self.SMILES))
            self.trans = ''

    def translate(self):
        form = self.br.select_form("form")
        form.set("smiles", self.SMILES)
        self.resp = self.br.submit_selected()
        soup = BeautifulSoup(self.resp.text, 'html.parser')
        if soup.title.text == 'USMILES Result':
            if soup.find('b') is not None:
                self.result = soup.find('b').text
                self.br.close()
                return self.result
            else:
                self.br.close()
                raise ValueError("[Parse] Cannot find results in the page returned.")
        elif soup.title.text == 'Translation Error':
            if soup.find('h3') is not None:
                self.br.close()
                raise ValueError(soup.find('h3').text)
            else:
                self.br.close()
                raise ValueError('[Translation] Unknown translation error.')
        else:
            self.br.close()
            raise ValueError('Unknown error.')