## Moments for the Angel sharks ##
{```{bash}```
mkdir WEA_sym WEA_asym AEW_sym AEW_asym null

python3 ~/bin/easySFS.py -i SNP.TRS.F07_DADI.recode.vcf -p AEW_popmap --preview
python3 ~/bin/easySFS.py -f -i SNP.TRS.F07_DADI.recode.vcf -o output_AEW -p AEW_popmap --proj 32,38,42
python3 ~/bin/easySFS.py -f -i SNP.TRS.F07_DADI.recode.vcf -o output_WEA -p WEA_popmap --proj 42,38,32

grep dDocent SNP.TRS.F07_DADI.recode.vcf | grep -v "#" | cut -f1 | sed 's/-/_/g' > loci.txt
sed -i 's/N//g' reference.fasta
python3 ~/bin/seq_length.py reference.fasta > ref_length.txt
grep -wf loci.txt ref_length.txt | cut -f2 | paste -sd+ | bc
} #notepad cleanup

#Plotting the models
{```{python}```
import sys
import os
import numpy
import scipy
import moments
sys.path.append('.')
import Angel_models
import Config_Make_fs_from_vcf as cfg

FILE = "output_AEW/dadi/Atlantic-East-West.sfs"
POPMAP = "AEW_popmap"
POP = ['Atl','East','West']

dd = moments.Misc.make_data_dict_vcf(FILE, POPMAP)
fs = moments.Spectrum.from_file(FILE)
ns = fs.sample_sizes

#print some useful information about the afs or jsfs
print("\n\n============================================================================")
print("\nData for site frequency spectrum\n")
print("Projection: {}".format(fs.pop_ids))
print("Sample sizes: {}".format(fs.sample_sizes))
print("Sum of SFS: {}".format(numpy.around(fs.S(), 2)))
print("\n============================================================================\n")

#Model_AEW_sym
model_func = Angel_models.pop_model_sym
params=[1,2,1,2,1,2,1,2,1,1,1,1,1]
model = moments.ModelPlot.generate_model(model_func, params, [10,10,10])

moments.ModelPlot.plot_model(model,
                             save_file='Model_AEW_sym.png',
                             fig_title='Example model from AEW with symetrical migration',
                             draw_scale=False,
                             pop_labels=fs.pop_ids,
                             nref=None,
                             gen_time=1.0,
                             gen_time_units='generations',
                             reverse_timeline=True)

#Model_AEW_asym
model_func = Angel_models.pop_model_asym
params=[1,2,1,2,1,2,1,2,1,0.5,1,0.5,1,0.5,1,1]
model = moments.ModelPlot.generate_model(model_func, params, [10,10,10])

moments.ModelPlot.plot_model(model,
                             save_file='Model_AEW_asym.png',
                             fig_title='Example model from AEW with asymetrical migration',
                             draw_scale=False,
                             pop_labels=fs.pop_ids,
                             nref=None,
                             gen_time=1.0,
                             gen_time_units='generations',
                             reverse_timeline=True)

#Model_null
model_func = Angel_models.null_model
params=[1,1,1,1]
model = moments.ModelPlot.generate_model(model_func, params, [10,10,10])

moments.ModelPlot.plot_model(model,
                             save_file='Model_null.png',
                             fig_title='Example of a null model with no migration',
                             draw_scale=False,
                             pop_labels=fs.pop_ids,
                             nref=None,
                             gen_time=1.0,
                             gen_time_units='generations',
                             reverse_timeline=True)

#Model_AEW_secondary contact and symetric migration
model_func = Angel_models.pop_model_second_sym
params=[1,2,1,2,1,2,1,2,1,1,1,1,1,1]
model = moments.ModelPlot.generate_model(model_func, params, [10,10,10])

moments.ModelPlot.plot_model(model,
                             save_file='Model_AEW_sec.png',
                             fig_title='Example model from AEW with symetrical migration and secondary contact',
                             draw_scale=False,
                             pop_labels=fs.pop_ids,
                             nref=None,
                             gen_time=1.0,
                             gen_time_units='generations',
                             reverse_timeline=True)
} #notepad cleanup

#### Workstatioms ####
## AEW Symetrical ##
{```{bash}```
#Round 1
seq 1 60 | sed 's/^/nudge/g' | xargs -I{} -P15 timeout 2d python3 Model_sym_Preturb.py {}
TMP=$(tail -n1 nudge*optimized.txt | sort -u | grep Round | wc -l); echo "$TMP iterations of 60 completed"

cat *optimized.txt | sort -k3,3n | uniq | sort -k1,1 -k3,3n > model_AEW.optimized.txt
tail -n1 model_AEW.optimized.txt > Best_run.txt
cp Best_run.txt Best_run_R1.txt
ls nudge*optimized.txt | while read i; do head -n1 $i >> $i; done

#Round 2
PARAMS=$(cut -f7 Best_run.txt)
sed -i "19s/.*/params=[$PARAMS]/" Model_sym_Preturb.py
sed -i 's/maxiters=\[3\]/maxiters=\[5\]/g' Config_Make_fs_from_vcf.py
sed -i 's/folds=\[3\]/folds=\[2\]/g' Config_Make_fs_from_vcf.py

seq 1 60 | sed 's/^/nudge/g' | xargs -I{} -P15 timeout 2d python3 Model_sym_Preturb.py {}
TMP=$(tail -n1 nudge*optimized.txt | sort -u | grep Round | wc -l); echo "$TMP iterations of 60 completed"

cat *optimized.txt | sort -k3,3n | uniq | sort -k1,1 -k3,3n > model_AEW.optimized.txt
tail -n1 model_AEW.optimized.txt > Best_run.txt
cp Best_run.txt Best_run_R2.txt
ls nudge*optimized.txt | while read i; do head -n1 $i >> $i; done

#Round 3
PARAMS=$(cut -f7 Best_run.txt)
sed -i "19s/.*/params=[$PARAMS]/" Model_sym_Preturb.py
sed -i 's/maxiters=\[5\]/maxiters=\[10\]/g' Config_Make_fs_from_vcf.py

seq 1 60 | sed 's/^/nudge/g' | xargs -I{} -P15 timeout 2d python3 Model_sym_Preturb.py {}
TMP=$(tail -n1 nudge*optimized.txt | sort -u | grep Round | wc -l); echo "$TMP iterations of 60 completed"

cat *optimized.txt | sort -k3,3n | uniq | sort -k1,1 -k3,3n > model_AEW.optimized.txt
tail -n1 model_AEW.optimized.txt > Best_run.txt
cp Best_run.txt Best_run_R3.txt
ls nudge*optimized.txt | while read i; do head -n1 $i >> $i; done

#Round 4
PARAMS=$(cut -f7 Best_run.txt)
sed -i "19s/.*/params=[$PARAMS]/" Model_sym_Preturb.py
sed -i 's/maxiters=\[10\]/maxiters=\[15\]/g' Config_Make_fs_from_vcf.py
sed -i 's/folds=\[2\]/folds=\[1\]/g' Config_Make_fs_from_vcf.py

seq 1 60 | sed 's/^/nudge/g' | xargs -I{} -P15 timeout 2d python3 Model_sym_Preturb.py {}
TMP=$(tail -n1 nudge*optimized.txt | sort -u | grep Round | wc -l); echo "$TMP iterations of 60 completed"

#Gathering output
cat *optimized.txt | sort -k3,3n | uniq | sort -k1,1 -k3,3n > model_AEW.optimized.txt
tail -n1 model_AEW.optimized.txt > Best_run.txt
cp Best_run.txt Best_run_R4.txt
ls nudge*optimized.txt | while read i; do head -n1 $i >> $i; done

#Resetting Config file
sed -i 's/maxiters=\[15\]/maxiters=\[3\]/g' Config_Make_fs_from_vcf.py
sed -i 's/folds=\[1\]/folds=\[3\]/g' Config_Make_fs_from_vcf.py
PARAMS=$(cut -f7 Best_run.txt)
sed -i "19s/.*/params=[$PARAMS]/" Model_sym_Preturb.py

#Killing current processes when necessary
for i in $(ps -fu afields | awk '/Model_sym_Preturb.py/ {print $2}'); do kill $i; done

#Restarting run where it left off
if [ -s rerun.txt ]; then rm rerun.txt; fi
seq 1 60 | while read i; do
if [ $(tail -n1 nudge$i.*optimized.txt | grep Round_1_Replicate_1 | wc -l) -gt 0 ]; then continue; fi
echo $i >> rerun.txt
done

cat rerun.txt | sed 's/^/nudge/g' | xargs -I{} -P15 timeout 2d python3 Model_sym_Preturb.py {}
} #notepad cleanup

## AEW Asymetrical ##
{```{bash}```
#Round 1
seq 1 60 | sed 's/^/nudge/g' | xargs -I{} -P15 timeout 2d python3 Model_asym_Preturb.py {}
TMP=$(tail -n1 nudge*optimized.txt | sort -u | grep Round | wc -l); echo "$TMP iterations of 60 completed"

cat *optimized.txt | sort -k3,3n | uniq | sort -k1,1 -k3,3n > model_AEW.optimized.txt
tail -n1 model_AEW.optimized.txt > Best_run.txt
cp Best_run.txt Best_run_R1.txt
ls nudge*optimized.txt | while read i; do head -n1 $i >> $i; done

#Round 2
PARAMS=$(cut -f7 Best_run.txt)
sed -i "19s/.*/params=[$PARAMS]/" Model_asym_Preturb.py
sed -i 's/maxiters=\[3\]/maxiters=\[5\]/g' Config_Make_fs_from_vcf.py
sed -i 's/folds=\[3\]/folds=\[2\]/g' Config_Make_fs_from_vcf.py

seq 1 60 | sed 's/^/nudge/g' | xargs -I{} -P15 timeout 2d python3 Model_asym_Preturb.py {}
TMP=$(tail -n1 nudge*optimized.txt | sort -u | grep Round | wc -l); echo "$TMP iterations of 60 completed"

cat *optimized.txt | sort -k3,3n | uniq | sort -k1,1 -k3,3n > model_AEW.optimized.txt
tail -n1 model_AEW.optimized.txt > Best_run.txt
cp Best_run.txt Best_run_R2.txt
ls nudge*optimized.txt | while read i; do head -n1 $i >> $i; done

#Round 3
PARAMS=$(cut -f7 Best_run.txt)
sed -i "19s/.*/params=[$PARAMS]/" Model_asym_Preturb.py
sed -i 's/maxiters=\[5\]/maxiters=\[10\]/g' Config_Make_fs_from_vcf.py

seq 1 60 | sed 's/^/nudge/g' | xargs -I{} -P15 timeout 2d python3 Model_asym_Preturb.py {}
TMP=$(tail -n1 nudge*optimized.txt | sort -u | grep Round | wc -l); echo "$TMP iterations of 60 completed"

cat *optimized.txt | sort -k3,3n | uniq | sort -k1,1 -k3,3n > model_AEW.optimized.txt
tail -n1 model_AEW.optimized.txt > Best_run.txt
cp Best_run.txt Best_run_R3.txt
ls nudge*optimized.txt | while read i; do head -n1 $i >> $i; done

#Round 4
PARAMS=$(cut -f7 Best_run.txt)
sed -i "19s/.*/params=[$PARAMS]/" Model_asym_Preturb.py
sed -i 's/maxiters=\[10\]/maxiters=\[15\]/g' Config_Make_fs_from_vcf.py
sed -i 's/folds=\[2\]/folds=\[1\]/g' Config_Make_fs_from_vcf.py

seq 1 60 | sed 's/^/nudge/g' | xargs -I{} -P15 timeout 2d python3 Model_asym_Preturb.py {}
TMP=$(tail -n1 nudge*optimized.txt | sort -u | grep Round | wc -l); echo "$TMP iterations of 60 completed"

#Gathering output
cat *optimized.txt | sort -k3,3n | uniq | sort -k1,1 -k3,3n > model_AEW.optimized.txt
tail -n1 model_AEW.optimized.txt > Best_run.txt
cp Best_run.txt Best_run_R4.txt
ls nudge*optimized.txt | while read i; do head -n1 $i >> $i; done

#Resetting Config file
sed -i 's/maxiters=\[15\]/maxiters=\[3\]/g' Config_Make_fs_from_vcf.py
sed -i 's/folds=\[1\]/folds=\[3\]/g' Config_Make_fs_from_vcf.py
PARAMS=$(cut -f7 Best_run.txt)
sed -i "19s/.*/params=[$PARAMS]/" Model_asym_Preturb.py

#Killing processes associated with a moments run
for i in $(ps -fu afields | awk '/Model_asym_Preturb.py/ {print $2}'); do kill $i; done

#Restarting run where it left off
if [ -s rerun.txt ]; then rm rerun.txt; fi
seq 1 60 | while read i; do
if [ $(tail -n1 nudge$i.*optimized.txt | grep Round_1_Replicate_1 | wc -l) -gt 0 ]; then continue; fi
echo $i >> rerun.txt
done

cat rerun.txt | sed 's/^/nudge/g' | xargs -I{} -P15 timeout 2d python3 Model_asym_Preturb.py {}
}

## WEA Asymetrical ##
{```{bash}```
#Round 1
seq 1 60 | sed 's/^/nudge/g' | xargs -I{} -P15 timeout 2d python3 Model_asym_Preturb.py {}
TMP=$(tail -n1 nudge*optimized.txt | sort -u | grep Round | wc -l); echo "$TMP iterations of 60 completed"

cat *optimized.txt | sort -k3,3n | uniq | sort -k1,1 -k3,3n > model_WEA.optimized.txt
tail -n1 model_WEA.optimized.txt > Best_run.txt
cp Best_run.txt Best_run_R1.txt
ls nudge*optimized.txt | while read i; do head -n1 $i >> $i; done

#Round 2
PARAMS=$(cut -f7 Best_run.txt)
sed -i "19s/.*/params=[$PARAMS]/" Model_asym_Preturb.py
sed -i 's/maxiters=\[3\]/maxiters=\[5\]/g' Config_Make_fs_from_vcf.py
sed -i 's/folds=\[3\]/folds=\[2\]/g' Config_Make_fs_from_vcf.py

seq 1 60 | sed 's/^/nudge/g' | xargs -I{} -P15 timeout 2d python3 Model_asym_Preturb.py {}
TMP=$(tail -n1 nudge*optimized.txt | sort -u | grep Round | wc -l); echo "$TMP iterations of 60 completed"

cat *optimized.txt | sort -k3,3n | uniq | sort -k1,1 -k3,3n > model_WEA.optimized.txt
tail -n1 model_WEA.optimized.txt > Best_run.txt
cp Best_run.txt Best_run_R2.txt
ls nudge*optimized.txt | while read i; do head -n1 $i >> $i; done

#Round 3
PARAMS=$(cut -f7 Best_run.txt)
sed -i "19s/.*/params=[$PARAMS]/" Model_asym_Preturb.py
sed -i 's/maxiters=\[5\]/maxiters=\[10\]/g' Config_Make_fs_from_vcf.py

seq 1 60 | sed 's/^/nudge/g' | xargs -I{} -P15 timeout 2d python3 Model_asym_Preturb.py {}
TMP=$(tail -n1 nudge*optimized.txt | sort -u | grep Round | wc -l); echo "$TMP iterations of 60 completed"

cat *optimized.txt | sort -k3,3n | uniq | sort -k1,1 -k3,3n > model_WEA.optimized.txt
tail -n1 model_WEA.optimized.txt > Best_run.txt
cp Best_run.txt Best_run_R3.txt
ls nudge*optimized.txt | while read i; do head -n1 $i >> $i; done

#Round 4
PARAMS=$(cut -f7 Best_run.txt)
sed -i "19s/.*/params=[$PARAMS]/" Model_asym_Preturb.py
sed -i 's/maxiters=\[10\]/maxiters=\[15\]/g' Config_Make_fs_from_vcf.py
sed -i 's/folds=\[2\]/folds=\[1\]/g' Config_Make_fs_from_vcf.py

seq 1 60 | sed 's/^/nudge/g' | xargs -I{} -P15 timeout 2d python3 Model_asym_Preturb.py {}
TMP=$(tail -n1 nudge*optimized.txt | sort -u | grep Round | wc -l); echo "$TMP iterations of 60 completed"

#Gathering output
cat *optimized.txt | sort -k3,3n | uniq | sort -k1,1 -k3,3n > model_WEA.optimized.txt
tail -n1 model_WEA.optimized.txt > Best_run.txt
cp Best_run.txt Best_run_R4.txt
ls nudge*optimized.txt | while read i; do head -n1 $i >> $i; done

#Killing current processes
for i in $(ps -fu afields | awk '/Model_asym_Preturb.py/ {print $2}'); do kill $i; done

#Resetting Config file
sed -i 's/maxiters=\[15\]/maxiters=\[3\]/g' Config_Make_fs_from_vcf.py
sed -i 's/folds=\[1\]/folds=\[3\]/g' Config_Make_fs_from_vcf.py
PARAMS=$(cut -f7 Best_run.txt)
sed -i "19s/.*/params=[$PARAMS]/" Model_asym_Preturb.py

#Restarting run where it left off
if [ -s rerun.txt ]; then rm rerun.txt; fi
seq 1 60 | while read i; do
if [ $(tail -n1 nudge$i.*optimized.txt | grep Round_1_Replicate_1 | wc -l) -gt 0 ]; then continue; fi
echo $i >> rerun.txt
done

cat rerun.txt | sed 's/^/nudge/g' | xargs -I{} -P15 timeout 2d python3 Model_asym_Preturb.py {}
} #notepad cleanup

## WEA Symetrical ##
{```{bash}```
#Round 1
seq 1 60 | sed 's/^/nudge/g' | xargs -I{} -P15 timeout 2d python3 Model_sym_Preturb.py {}
TMP=$(tail -n1 nudge*optimized.txt | sort -u | grep Round | wc -l); echo "$TMP iterations of 60 completed"

cat *optimized.txt | sort -k3,3n | uniq | sort -k1,1 -k3,3n > model_WEA.optimized.txt
tail -n1 model_WEA.optimized.txt > Best_run.txt
cp Best_run.txt Best_run_R1.txt
ls nudge*optimized.txt | while read i; do head -n1 $i >> $i; done

#Round 2
PARAMS=$(cut -f7 Best_run.txt)
sed -i "19s/.*/params=[$PARAMS]/" Model_sym_Preturb.py
sed -i 's/maxiters=\[3\]/maxiters=\[5\]/g' Config_Make_fs_from_vcf.py
sed -i 's/folds=\[3\]/folds=\[2\]/g' Config_Make_fs_from_vcf.py

seq 1 60 | sed 's/^/nudge/g' | xargs -I{} -P15 timeout 2d python3 Model_sym_Preturb.py {}
TMP=$(tail -n1 nudge*optimized.txt | sort -u | grep Round | wc -l); echo "$TMP iterations of 60 completed"

cat *optimized.txt | sort -k3,3n | uniq | sort -k1,1 -k3,3n > model_WEA.optimized.txt
tail -n1 model_WEA.optimized.txt > Best_run.txt
cp Best_run.txt Best_run_R2.txt
ls nudge*optimized.txt | while read i; do head -n1 $i >> $i; done

#Round 3
PARAMS=$(cut -f7 Best_run.txt)
sed -i "19s/.*/params=[$PARAMS]/" Model_sym_Preturb.py
sed -i 's/maxiters=\[5\]/maxiters=\[10\]/g' Config_Make_fs_from_vcf.py

seq 1 60 | sed 's/^/nudge/g' | xargs -I{} -P15 timeout 2d python3 Model_sym_Preturb.py {}
TMP=$(tail -n1 nudge*optimized.txt | sort -u | grep Round | wc -l); echo "$TMP iterations of 60 completed"

cat *optimized.txt | sort -k3,3n | uniq | sort -k1,1 -k3,3n > model_WEA.optimized.txt
tail -n1 model_WEA.optimized.txt > Best_run.txt
cp Best_run.txt Best_run_R3.txt
ls nudge*optimized.txt | while read i; do head -n1 $i >> $i; done

#Round 4
PARAMS=$(cut -f7 Best_run.txt)
sed -i "19s/.*/params=[$PARAMS]/" Model_sym_Preturb.py
sed -i 's/maxiters=\[10\]/maxiters=\[15\]/g' Config_Make_fs_from_vcf.py
sed -i 's/folds=\[2\]/folds=\[1\]/g' Config_Make_fs_from_vcf.py

seq 1 60 | sed 's/^/nudge/g' | xargs -I{} -P15 timeout 2d python3 Model_sym_Preturb.py {}
TMP=$(tail -n1 nudge*optimized.txt | sort -u | grep Round | wc -l); echo "$TMP iterations of 60 completed"

cat *optimized.txt | sort -k3,3n | uniq | sort -k1,1 -k3,3n > model_WEA.optimized.txt
tail -n1 model_WEA.optimized.txt > Best_run.txt
cp Best_run.txt Best_run_R4.txt
ls nudge*optimized.txt | while read i; do head -n1 $i >> $i; done

#Killing current processes
for i in $(ps -fu afields | awk '/Model_sym_Preturb.py/ {print $2}'); do kill $i; done

#Resetting Config file
sed -i 's/maxiters=\[15\]/maxiters=\[3\]/g' Config_Make_fs_from_vcf.py
sed -i 's/folds=\[1\]/folds=\[3\]/g' Config_Make_fs_from_vcf.py
PARAMS=$(cut -f7 Best_run.txt)
sed -i "19s/.*/params=[$PARAMS]/" Model_asym_Preturb.py
} #notepad cleanup

## Null ##
{```{bash}```
#Round 1
seq 1 60 | sed 's/^/nudge/g' | xargs -I{} -P15 timeout 2d python3 Model_null_Preturb.py {}

cat *optimized.txt | sort -k3,3n | uniq | sort -k1,1 -k3,3n > model_null.optimized.txt
tail -n1 model_null.optimized.txt > Best_run.txt
cp Best_run.txt Best_run_R1.txt
ls nudge*optimized.txt | while read i; do head -n1 $i >> $i; done

#Round 2
PARAMS=$(cut -f7 Best_run.txt)
sed -i "19s/.*/params=[$PARAMS]/" Model_null_Preturb.py
sed -i 's/maxiters=\[3\]/maxiters=\[5\]/g' Config_Make_fs_from_vcf.py
sed -i 's/folds=\[3\]/folds=\[2\]/g' Config_Make_fs_from_vcf.py

seq 1 60 | sed 's/^/nudge/g' | xargs -I{} -P15 timeout 2d python3 Model_null_Preturb.py {}

cat *optimized.txt | sort -k3,3n | uniq | sort -k1,1 -k3,3n > model_null.optimized.txt
tail -n1 model_null.optimized.txt > Best_run.txt
cp Best_run.txt Best_run_R2.txt
ls nudge*optimized.txt | while read i; do head -n1 $i >> $i; done

#Round 3
PARAMS=$(cut -f7 Best_run.txt)
sed -i "19s/.*/params=[$PARAMS]/" Model_null_Preturb.py
sed -i 's/maxiters=\[5\]/maxiters=\[10\]/g' Config_Make_fs_from_vcf.py

seq 1 60 | sed 's/^/nudge/g' | xargs -I{} -P15 timeout 2d python3 Model_null_Preturb.py {}

cat *optimized.txt | sort -k3,3n | uniq | sort -k1,1 -k3,3n > model_null.optimized.txt
tail -n1 model_null.optimized.txt > Best_run.txt
cp Best_run.txt Best_run_R3.txt
ls nudge*optimized.txt | while read i; do head -n1 $i >> $i; done

#Round 4
PARAMS=$(cut -f7 Best_run.txt)
sed -i "19s/.*/params=[$PARAMS]/" Model_null_Preturb.py
sed -i 's/maxiters=\[10\]/maxiters=\[15\]/g' Config_Make_fs_from_vcf.py
sed -i 's/folds=\[2\]/folds=\[1\]/g' Config_Make_fs_from_vcf.py

seq 1 60 | sed 's/^/nudge/g' | xargs -I{} -P15 timeout 2d python3 Model_null_Preturb.py {}

cat *optimized.txt | sort -k3,3n | uniq | sort -k1,1 -k3,3n > model_null.optimized.txt
tail -n1 model_null.optimized.txt > Best_run.txt
cp Best_run.txt Best_run_R4.txt
ls nudge*optimized.txt | while read i; do head -n1 $i >> $i; done

#Resetting Config file
sed -i 's/maxiters=\[15\]/maxiters=\[3\]/g' Config_Make_fs_from_vcf.py
sed -i 's/folds=\[1\]/folds=\[3\]/g' Config_Make_fs_from_vcf.py
} #notepad cleanup

#Getting CI (AEW_sym) 
{```{python3}```
import sys
import os
import moments
import time
import hickle
sys.path.append('~/Workspace/Angels/moments2')
import Angel_models
sys.path.append('.')
import Config_Make_fs_from_vcf as cfg

model = Angel_models.pop_model_sym
params=[1.8632,191.2111,1.9303,1.2138,0.9394,6644.2825,14.6523,267.6575,0.0255,28.807,57.1529,9.8295,24.4424]
fs = moments.Spectrum.from_file(cfg.FILE)
ns = fs.sample_sizes
dd = moments.Misc.make_data_dict_vcf("../SNP.TRS.F07_DADI.recode.vcf", cfg.POPMAP)

import time
start = time.process_time()
func_ex = model
all_boot = moments.Misc.bootstrap(dd, fs.pop_ids, ns, mask_corners=True, polarized=False, num_boots=1000)
p0 = params
data = fs
hickle.dump(all_boot, 'all_boot_GIM_Feb22.hkl.gz', mode='w', compression='gzip')
print(time.process_time() - start)

start = time.process_time()
GIM =  moments.Godambe.GIM_uncert(func_ex, all_boot, p0, data, log=False, multinom=True, eps=0.01, return_GIM=True)
hickle.dump(GIM, 'GIM_Feb22.hkl.gz', mode='w', compression='gzip')
print(time.process_time() - start)

GIM[0]

start = time.process_time()
all_boot = moments.Misc.bootstrap(dd, fs.pop_ids, ns, mask_corners=True, polarized=False, num_boots=1000)
hickle.dump(all_boot, 'all_boot_GIM_May07.hkl.gz', mode='w', compression='gzip')
print(time.process_time() - start)

start = time.process_time()
GIM =  moments.Godambe.GIM_uncert(func_ex, all_boot, p0, data, log=False, multinom=True, eps=0.01, return_GIM=True)
hickle.dump(GIM, 'GIM_Mar07.hkl.gz', mode='w', compression='gzip')
print(time.process_time() - start)

GIM[0]
} #notepad cleanup

#Getting AIC values for a given set of parameters and model
{```{python}```
import sys
import os
import numpy
import scipy
import moments
sys.path.append('..')
import Angel_models
sys.path.append('.')
import Config_Make_fs_from_vcf as cfg

FILE = "../output_WEA/dadi/West-East-Atlantic.sfs"
POPMAP = "WEA_popmap"
POP = ['Atl','East','West']

POP = ['West','East','Atl']

dd = moments.Misc.make_data_dict_vcf(FILE, POPMAP)
fs = moments.Spectrum.from_file(FILE)
ns = fs.sample_sizes

#print some useful information about the afs or jsfs
print("\n\n============================================================================")
print("\nData for site frequency spectrum\n")
print("Projection: {}".format(fs.pop_ids))
print("Sample sizes: {}".format(fs.sample_sizes))
print("Sum of SFS: {}".format(numpy.around(fs.S(), 2)))
print("\n============================================================================\n")

#Model_AEW_sym
model_func = Angel_models.pop_model_sym
params=[13.4626, 11.9305, 16.4862, 199.4982, 0.0223, 0.3551, 0.1347, 0.9797, 97.9553, 96.6461, 16.1688, 14.3861, 0.3624]

model = moments.ModelPlot.generate_model(model_func, params, [10,10,10])

model_fit = model_func(params, fs.sample_sizes)
ll = moments.Inference.ll(model_fit, fs)
AIC_value = ( -2*( float(ll))) + (2*len(params))
} #notepad cleanup
