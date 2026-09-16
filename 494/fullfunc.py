import structscrape
import layer1_seqanalysis
import layer2_structgeometry
import layer3_electrostatics
import ast

'''
wanna run all in cmd with inputs, start by prompting the user to input uniprot search terms, operator, identity

a) run uniprot search command from structscrape, download all structures to a working folder, the rest of data should be on excel
including references for pdb structures

b) run created excel sheet through all 3 layers in the pipeline, paste all ideal traits into one big dataframe

c) pull all that shit, perform pca and output
'''


def run_tool():

    sterms = ast.literal_eval(input("Please input a python list of search terms (maybe i need to do a dict): "))
    ops = input('Please input your operator for your search term (AND or OR): ')
    id = float(input("Input identity for redundancy filter: "))

    metadata = structscrape.uniprot_search(sterms, id, ops)