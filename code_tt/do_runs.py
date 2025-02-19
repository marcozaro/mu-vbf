
import run as run
import time
import cluster as cluster
import os

me_dir = os.getcwd()

# parse a reasonable run_card
input_dict = run.parse_input('run_card_checkmassive_ymax.dat')

input_dict['ecm'] = 10000
runs = []
print('STARTING RUNS')


for ecm in [3000, ]:
    for muf in [10]:
        for ymin in [-5, -2.5]:

            input_dict['fixscale'] = muf
            input_dict['ecm'] = ecm
            input_dict['ymin'] = ymin

            run_name = 'run_TEST_%dtev_comparemassive_ymax%3.1f_mufix%d' % (ecm,-ymin,muf)

            runs.append(run_name)
            print(input_dict, run_name)
            run.write_input('input.inc', input_dict)
            run.write_printout('printout.inc', input_dict)
            run.compile()
            run.run(run_name, input_dict, prepareonly=True)


print(runs)


#import multiprocessing as mp
#cpus = mp.cpu_count()
#pool = mp.Pool(processes=cpus)

run_cluster = cluster.MultiCore(8)

for rr in runs:
    run_cluster.submit2("./driver", [], cwd = rr)
                            

update_status = lambda i, r, f: print('IDLE %d; RUNNING %d; DONE %d' % (i,r,f))
try:
    run_cluster.wait(me_dir, update_status)
except:
    run_cluster.remove()
    raise



#rr = pool.map(run.run_folder, runs)
