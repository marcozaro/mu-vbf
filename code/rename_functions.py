#!/usr/bin/python

import sys

infile = sys.argv[1]
suffix = sys.argv[2]
content = open(infile).read()

outfile = open(infile.replace('.f', '_mod.f'), 'w')

# read the subroutines and functions
subroutines = []
functions = []
for l in content.split('\n'):
    if not l.strip() or l.strip()[0] in ['C', '!']:
        continue
    if l.strip().startswith('END'):
        continue
    if l.strip().startswith('SUBROUTINE'):
        subroutines.append(l.split('SUBROUTINE')[1].split('(')[0].strip())
    elif "FUNCTION" in l:
        functions.append(l.split('FUNCTION')[1].split('(')[0].strip())



for s in subroutines:
    content = content.replace('SUBROUTINE %s(' % s, 'SUBROUTINE %s_%s(' %(s, suffix))
    content = content.replace('CALL %s(' % s, 'CALL %s_%s(' %(s, suffix))
    content = content.replace('END SUBROUTINE %s\n' % s, 'END SUBROUTINE\n')

for f in functions:
    #content = content.replace('FUNCTION %s(' % f, 'FUNCTION %s_%s(' %(f, suffix))
    content = content.replace('%s' % f, '%s_%s' %(f, suffix))


common_blocks = ['ALL_SQSPLITORDERS', 
                 'AMPSPLITORDERS',
                 'COLOR_CORRELATION_MAPS' ,
                 'COLOR_CORRELATIONS' ,
                 'BORN_HEL_CONFIGS' ,
                 'CHOSEN_BORN_SQSO' ,
                 'COLOR_CONNECTION_DEFINITIONS_NLO' ,
                 'CC_INDEX_TO_DEFINITION_NLO' ,
                 'BORN_BEAM_POL' ,
                 'HELUSERCHOICE' ,
                 'SPIN_CORRELATION_DATA' ,
                 'ALL_SQSPLITORDERS' ,
                 'AMPSPLITORDERS' ,
                 'COLOR_CORRELATION_MAPS' ,
                 'COLOR_CORRELATIONS' ,
                 'BORN_HEL_CONFIGS' ,
                 'CHOSEN_BORN_SQSO' ,
                 'COLOR_CONNECTION_DEFINITIONS_NLO' ,
                 'CC_INDEX_TO_DEFINITION_NLO' ,
                 'BORN_BEAM_POL' ,
                 'HELUSERCHOICE' ,
                 'SPIN_CORRELATION_DATA' ,
                 'ALL_SQSPLITORDERS' ,
                 'AMPSPLITORDERS' ,
                 'BORN_HEL_CONFIGS', 
                 'CHOSEN_BORN_SQSO' ,
                 'BORN_BEAM_POL' ,
                 'HELUSERCHOICE' ]

for c in common_blocks:
    content =  content.replace('COMMON/%s/' % c, 'COMMON/%s_%s/' %(c, suffix))

outfile.write(content)
outfile.close()






