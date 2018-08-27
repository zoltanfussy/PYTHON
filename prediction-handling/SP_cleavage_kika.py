handle = 'further'
inFile = open(handle + '.tsv', 'r')
#input file v tabulke nazov_kontigu, SP_cleavage, prot_seq

for protein in inFile:
##    rozdelenie riadku po jednotlivych bunkach v tabulke
    splitted_line = protein.split('\t')
    contig = splitted_line[0]
##    TMD = int(splitted_line[1])
    SP = splitted_line[1]
    seq = splitted_line[2]
##v proteinoch obsahujucich SP obstiepi N-koniec podla dlzky SP
    try:
        SP = int(SP)
        SP_cleaved = '>{}\n{}'.format(contig,seq[SP:])
##        print(SP_cleaved)
        with open(handle + '_cleaved.txt', 'a') as result:
            result.write(SP_cleaved)
##contig s proteinom, ktory neobsahuje SP, zapise do suboru s chybami
    except ValueError as VE:
        with open(handle + 'ValueErrors.txt', 'a') as ValueErrors:
            ValueErrors.write('{}\t{}\n'.format(contig,str(VE)))
