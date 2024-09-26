
import run as run
import time

# parse a reasonable run_card
input_dict = run.parse_input('run_card_1.dat')

###MZMZ no mll cut
input_dict['mllmin'] = 0.

for ecm in [1000, 3000, 10000]:
    for ymin in [-2.5, -5]:
        input_dict['ecm'] = ecm
        input_dict['ymin'] = ymin

        run_name = 'run_nlo_%d_nomcut_ymin%3.1f' % ( ecm/1000., ymin)

        print(input_dict, run_name)
        run.write_input('input.inc', input_dict)
        run.write_printout('printout.inc', input_dict)
        run.compile()
        run.run(run_name, input_dict)
